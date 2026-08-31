/* =====================================================================
   script.js -- form handling, output rendering and the self-test.
   All the mathematics lives in nerf.js, which this file assumes is
   already loaded.
   ===================================================================== */
"use strict";

const $ = (id) => document.getElementById(id);

/* =====================================================================
   Small helpers
   ===================================================================== */

/** Read a numeric field, or throw a message naming the field. */
function num(id, label) {
  const v = parseFloat($(id).value);
  if (!isFinite(v)) throw label + ' is not a number ("' + $(id).value + '")';
  return v;
}

/** Show or clear a Bootstrap alert. */
function setError(id, message) {
  const el = $(id);
  el.textContent = message || "";
  el.classList.toggle("d-none", !message);
}

/** The four computed positions, as a table. Atom 4 is emphasised. */
function xyzTable(P) {
  const rows = P.map(
    (p, i) =>
      "<tr" +
      (i === 3 ? ' class="fw-bold"' : "") +
      "><td>" +
      (i + 1) +
      "</td>" +
      p.map((v) => '<td class="text-end">' + v.toFixed(4) + "</td>").join("") +
      "</tr>",
  ).join("");

  return (
    '<table class="table table-sm font-monospace small tabular mt-3 mb-0">' +
    '<thead><tr class="text-secondary">' +
    '<th scope="col">atom</th>' +
    '<th scope="col" class="text-end">x (&Aring;)</th>' +
    '<th scope="col" class="text-end">y (&Aring;)</th>' +
    '<th scope="col" class="text-end">z (&Aring;)</th>' +
    "</tr></thead><tbody>" +
    rows +
    "</tbody></table>"
  );
}

/** Atom 4's internal coordinates measured back out of the computed XYZ. */
function checkBlock(P, r, theta, phi) {
  const line = (name, got, asked) =>
    "  " +
    name.padEnd(6) +
    "= " +
    got.toFixed(6).padStart(11) +
    "   (asked for " +
    asked.toFixed(6) +
    ")";

  return (
    '<pre class="bg-body-tertiary border rounded p-3 small mt-3 mb-0">' +
    "check, measured back from the result:\n" +
    line("r34", distOf(P[3], P[2]), r) +
    "\n" +
    line("theta", angOf(P[1], P[2], P[3]), theta) +
    "\n" +
    line("phi", dihOf(P[0], P[1], P[2], P[3]), wrap180(phi)) +
    "</pre>"
  );
}

/* =====================================================================
   Build all four atoms from internal coordinates
   ===================================================================== */
function runBuild() {
  try {
    const r12 = num("r12", "r12");
    const r23 = num("r23", "r23");
    const a123 = num("a123", "theta123");
    const r34 = num("r34", "r34");
    const a234 = num("a234", "theta234");
    const phi = num("d1234", "phi1234");

    if (r12 <= 0 || r23 <= 0 || r34 <= 0) throw "bond lengths must be positive";

    const P = buildFour(r12, r23, a123, r34, a234, phi);

    setError("berr", "");
    $("bout").innerHTML = xyzTable(P) + checkBlock(P, r34, a234, phi);
  } catch (e) {
    setError("berr", String(e));
    $("bout").innerHTML = "";
  }
}

/* =====================================================================
   Wiring
   ===================================================================== */
$("anti").addEventListener("click", () => {
  $("d1234").value = "180.0";
  runBuild();
});
$("gauche").addEventListener("click", () => {
  $("d1234").value = "60.0";
  runBuild();
});

["r12", "r23", "a123", "r34", "a234", "d1234"].forEach((id) =>
  $(id).addEventListener("input", runBuild),
);

/* =====================================================================
   Self-test -- runs on load, result logged to the console.
   ===================================================================== */
function selfTest() {
  const results = [];
  const check = (name, fn) => {
    try {
      results.push([name, fn()]);
    } catch (e) {
      results.push([name, String(e)]);
    }
  };

  // 1. Round trip. Build from internal coordinates, measure them back out of the
  //    Cartesian result, and require a match. This is what pins down the torsion
  //    sign convention between placeAtom() and dihOf().
  check("round trip", () => {
    const cases = [
      [1.526, 1.526, 112.7, 1.526, 112.7, 180],
      [1.09, 1.43, 109.5, 0.97, 104.5, -60],
      [1.33, 1.5, 120.0, 1.1, 118.0, 0],
      [2.05, 1.21, 175.0, 1.44, 91.3, 137.8],
      [1.0, 1.0, 90.0, 1.0, 90.0, -179.9],
    ];
    let worst = 0;

    for (const c of cases) {
      const P = buildFour.apply(null, c);
      worst = Math.max(
        worst,
        Math.abs(distOf(P[1], P[0]) - c[0]),
        Math.abs(distOf(P[2], P[1]) - c[1]),
        Math.abs(angOf(P[0], P[1], P[2]) - c[2]),
        Math.abs(distOf(P[3], P[2]) - c[3]),
        Math.abs(angOf(P[1], P[2], P[3]) - c[4]),
        angDiff(dihOf(P[0], P[1], P[2], P[3]), c[5]),
      );
    }
    return worst < 1e-9 ? true : "max deviation " + worst.toExponential(2);
  });

  // 2. Known geometry. Four bonds at the exact tetrahedral angle arccos(-1/3)
  //    must give six identical outer distances, each equal to bond * sqrt(8/3).
  check("tetrahedral geometry", () => {
    const r = 1.087,
      tet = Math.acos(-1 / 3) * R2D,
      C = [0, 0, 0];

    const H = [[r, 0, 0]];
    H.push(placeThird(H[0], C, r, tet));
    H.push(placeAtom(H[1], H[0], C, r, tet, 120));
    H.push(placeAtom(H[1], H[0], C, r, tet, -120));

    const d = [];
    for (let i = 0; i < 4; i++)
      for (let j = i + 1; j < 4; j++) d.push(distOf(H[i], H[j]));

    const spread = Math.max.apply(null, d) - Math.min.apply(null, d);
    const expected = r * Math.sqrt(8 / 3);

    if (spread > 1e-12) return "outer spread " + spread.toExponential(2);
    if (Math.abs(d[0] - expected) > 1e-12)
      return "got " + d[0].toFixed(9) + ", expected " + expected.toFixed(9);
    return true;
  });

  // 3. Rigid-body invariance. Rotating and translating atoms 1-3 must carry
  //    atom 4 with them: internal coordinates cannot depend on the frame the
  //    molecule happens to sit in.
  check("rigid-body invariance", () => {
    const A = [0.11, -0.42, 0.77],
      B = [1.63, 0.05, -0.21],
      C = [2.02, 1.44, 0.63];
    const D = placeAtom(A, B, C, 1.526, 112.7, -73.4);

    const axis = [0.37, -0.81, 0.45],
      ang = 1.234,
      shift = [-3.1, 7.7, 0.9];
    const T = (p) => add(rodrigues(p, axis, ang), shift);

    const moved = placeAtom(T(A), T(B), T(C), 1.526, 112.7, -73.4);
    const e = distOf(moved, T(D));
    return e < 1e-12 ? true : "moved by " + e.toExponential(2);
  });

  // 4. Atom 3 matches the closed form p3 = (r12 - r cos0, r sin0, 0).
  check("atom 3 closed form", () => {
    const r12 = 1.526,
      r = 1.43,
      th = 109.5;
    const p = placeThird([0, 0, 0], [r12, 0, 0], r, th);
    const want = [r12 - r * Math.cos(th * D2R), r * Math.sin(th * D2R), 0];

    const e = distOf(p, want);
    return e < 1e-12 ? true : "differs by " + e.toExponential(2);
  });

  const passed = results.filter((x) => x[1] === true).length;

  const report = results
    .map((x) =>
      x[1] === true ? "PASS  " + x[0] : "FAIL  " + x[0] + " -- " + x[1],
    )
    .join("\n");

  console.log(
    "self-test " + passed + "/" + results.length + "\n" + report,
  );

  return results;
}

/* boot */
runBuild();
selfTest();
