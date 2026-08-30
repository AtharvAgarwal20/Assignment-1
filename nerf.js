/* =====================================================================
   nerf.js -- the mathematical model, with no reference to the DOM.
   Internal coordinates (bond length, bond angle, torsion angle) to
   Cartesian coordinates, by the Natural Extension Reference Frame
   construction. The equations are set out in formulae.pdf.

   Plain script rather than an ES module: ES modules are blocked by CORS
   when a page is opened straight from the filesystem, and this has to
   work by double-clicking index.html.
   ===================================================================== */
"use strict";

/* ---------------------------------------------------------------------
   Vector helpers. A vector is a plain array [x, y, z].
   --------------------------------------------------------------------- */
const sub   = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add   = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const scl   = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot   = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1],
                         a[2] * b[0] - a[0] * b[2],
                         a[0] * b[1] - a[1] * b[0]];
const len   = a => Math.hypot(a[0], a[1], a[2]);
const unit  = a => scl(a, 1 / len(a));

const D2R = Math.PI / 180;
const R2D = 180 / Math.PI;

/* Rotate v about the axis k by ang radians (Rodrigues' formula).
   Only used by the rigid-body self-test. */
const rodrigues = (v, k, ang) => {
  k = unit(k);
  const c = Math.cos(ang), s = Math.sin(ang);
  return add(add(scl(v, c), scl(cross(k, v), s)), scl(k, dot(k, v) * (1 - c)));
};


/* =====================================================================
   Forward: internal coordinates -> Cartesian
   ===================================================================== */

/**
 * The local orthonormal frame at atom c, built from the three reference atoms.
 *   e1  along the b -> c bond          (local +x)
 *   e3  normal to the a-b-c plane      (local +z, the axis the torsion turns about)
 *   e2  e3 x e1, completing the set    (local +y)
 *
 * These are the columns of the rotation matrix M in formulae.pdf. The same
 * frame is used by dihOf() below, which is what keeps the forward and inverse
 * torsion sign conventions consistent with each other.
 *
 * @returns {Array<number[]>} [e1, e2, e3]
 */
function frame(pa, pb, pc) {
  const e1 = unit(sub(pc, pb));
  let n = cross(sub(pb, pa), e1);

  if (len(n) < 1e-10) {
    // a, b and c are collinear, so the a-b-c plane -- and with it the torsion
    // -- is undefined. Any perpendicular to e1 gives valid geometry; it only
    // fixes the molecule's arbitrary overall orientation.
    n = cross(e1, Math.abs(e1[1]) < 0.9 ? [0, 1, 0] : [0, 0, 1]);
  }

  const e3 = unit(n);
  return [e1, cross(e3, e1), e3];
}

/**
 * Place an atom given all three internal coordinates:
 *
 *   p = pc + M . ( -r cos0,  r sin0 cos1,  r sin0 sin1 )
 *
 * @param {number[]} pa  torsion reference atom
 * @param {number[]} pb  angle reference atom
 * @param {number[]} pc  bonded reference atom -- the new atom hangs off this one
 * @param {number}   r         bond length to pc, in Angstrom
 * @param {number}   thetaDeg  bond angle b-c-new, in degrees
 * @param {number}   phiDeg    torsion a-b-c-new, in degrees
 * @returns {number[]} the new atom's Cartesian position
 */
function placeAtom(pa, pb, pc, r, thetaDeg, phiDeg) {
  const [e1, e2, e3] = frame(pa, pb, pc);
  const th = thetaDeg * D2R;
  const ph = phiDeg * D2R;

  return add(pc, add(scl(e1, -r * Math.cos(th)),
                 add(scl(e2,  r * Math.sin(th) * Math.cos(ph)),
                     scl(e3,  r * Math.sin(th) * Math.sin(ph)))));
}

/**
 * Place atom 3: bond length and bond angle, but no torsion yet. Passing pb as
 * its own torsion reference sends frame() down the collinear branch, which
 * supplies the arbitrary perpendicular that a torsion would otherwise fix.
 * Reduces to the closed form p3 = (r12 - r cos0, r sin0, 0) in the canonical
 * case where atom 1 is at the origin and atom 2 lies on +x.
 */
const placeThird = (pb, pc, r, thetaDeg) => placeAtom(pb, pb, pc, r, thetaDeg, 0);

/**
 * The whole four-atom conversion.
 * @returns {Array<number[]>} [p1, p2, p3, p4]
 */
function buildFour(r12, r23, a123, r34, a234, d1234) {
  const p1 = [0, 0, 0];                                 // no data given: the origin
  const p2 = [r12, 0, 0];                               // bond length only: lay it on +x
  const p3 = placeThird(p1, p2, r23, a123);             // + bond angle: choose the xy-plane
  const p4 = placeAtom(p1, p2, p3, r34, a234, d1234);   // + torsion: the general case

  return [p1, p2, p3, p4];
}


/* =====================================================================
   Inverse: Cartesian -> internal coordinates.
   Used to check the forward result, and by the self-test.
   ===================================================================== */

/** Distance between two atoms. */
const distOf = (p, q) => len(sub(p, q));

/** Bond angle at q, i.e. the angle p-q-s, in degrees. */
const angOf = (p, q, s) =>
  Math.acos(Math.min(1, Math.max(-1, dot(unit(sub(p, q)), unit(sub(s, q)))))) * R2D;

/**
 * Signed torsion angle a-b-c-i, in degrees.
 * Reusing frame() makes this an exact inverse of placeAtom(): the bracketed
 * vector there has components (r sin0 cos1) along e2 and (r sin0 sin1) along e3,
 * so atan2 of those two recovers the torsion with the right sign, every time.
 */
function dihOf(pa, pb, pc, pi) {
  const [, e2, e3] = frame(pa, pb, pc);
  const v = sub(pi, pc);
  return Math.atan2(dot(v, e3), dot(v, e2)) * R2D;
}

/** Smallest absolute difference between two angles, in degrees. */
const angDiff = (a, b) => Math.abs(((a - b) % 360 + 540) % 360 - 180);

/** Fold an angle into (-180, 180], matching what atan2 reports. */
const wrap180 = a => {
  const n = ((a % 360) + 540) % 360 - 180;
  return n === -180 ? 180 : n;
};
