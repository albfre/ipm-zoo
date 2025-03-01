let wasmModule;

// When the page loads, initialize the module
document.addEventListener('DOMContentLoaded', function() {
  // Initialize the WASM module with debug logging
  console.log("Initializing WASM module...");
  SymbolicOptimizationModule().then(module => {
    console.log("WASM module initialized:", module);
    // Log all properties and methods of the module
    console.log("Module properties:", Object.getOwnPropertyNames(module));
    wasmModule = module;
    // Now that the module is ready, we can update the problem display
    updateProblem();
  }).catch(error => {
    console.error("Failed to initialize WASM module:", error);
  });
  
  // Set up event listeners for inputs
  document.querySelectorAll('input').forEach(input => {
    input.addEventListener('change', updateProblem);
  });
});

const x = "x"; // Primary variable
const A_ineq = "A"; // Inequality constraint matrix
const l_A = "l_{A}"; // Inequality constraint lower bound
const u_A = "u_{A}"; // Inequality constraint upper bound
const A_eq = "A_{\\text{eq}}"; // Equality constraint matrix
const b_eq = "b"; // Equality constraint right-hand side
const l_x = "l_x"; // Primary variable lower bound
const u_x = "u_x"; // Primary variable upper bound
const s_A = "s"; // Slack variable for inequality constraint
const s_lA = "g"; // Slack variable for inequality constraint lower bound
const s_uA = "t"; // Slack variable for inequality constraint upper bound
const s_lx = "y"; // Slack variable for x lower bound
const s_ux = "z"; // Slack variable for x upper bound
const getObjective = (logBarrierVariables = []) => {
  let output = "\\underset{" + x + "}{\\text{minimize}} \\quad & \\frac{1}{2} " + x + "^T Q " + x + " + c^T" + x;
  for (const v of logBarrierVariables) {
    output += "- \\mu e^T \\log(" + v + ")";
  }
  output += " \\\\";
  return output;
}

function getOriginalProblem(inequalities, equalities, variableBounds) {
  const hasInequalities = inequalities !== "none";
  const hasEqualities = equalities;
  const hasVariableBounds = variableBounds !== "none";
  let outputText = "";
  outputText += "\\[ \\begin{align*}\n";
  outputText += getObjective();

  if (hasInequalities || hasEqualities || hasVariableBounds) {
    outputText += "\\text{subject to} \\quad";
  }

  if (hasEqualities) {
    outputText += "& " + A_eq + x + " = " + b_eq + "\\\\\n"
  }

  if (hasInequalities) {
    if (inequalities === "lower") outputText += "& " + A_ineq + x + " \\geq " + l_A + " \\\\\n";
    if (inequalities === "upper") outputText += "& " + A_ineq + x + " \\leq " + u_A + " \\\\\n";
    if (inequalities === "both") outputText += "& " + l_A + " \\leq " + A_ineq + x + " \\leq " + u_A + " \\\\\n";
  }

  if (variableBounds !== "none") {
    if (variableBounds === "lower") outputText += "& " + x + " \\geq " + l_x + " \\\\\n";
    if (variableBounds === "upper") outputText += "& " + x + " \\leq " + u_x + " \\\\\n";
    if (variableBounds === "both") outputText += "& " + l_x + " \\leq " + x + " \\leq " + u_x + " \\\\\n";
  }

  outputText += "\\end{align*} \\]\n"
  return outputText;
}

function getNonnegativeSlacks(inequalities, equalities, variableBounds) {
  let nonnegativeSlacks = [];
  if (inequalities === "lower" || inequalities == "both") {
    nonnegativeSlacks.push(s_lA);
  }
  if (inequalities === "upper" || inequalities == "both") {
    nonnegativeSlacks.push(s_uA);
  }
  if (variableBounds === "lower" || variableBounds == "both") {
    nonnegativeSlacks.push(s_lx);
  }
  if (variableBounds === "upper" || variableBounds == "both") {
    nonnegativeSlacks.push(s_ux);
  }
  return nonnegativeSlacks;
}

function getSlackProblem(inequalities, equalities, variableBounds, logBarriers) {
  const hasInequalities = inequalities !== "none";
  const hasEqualities = equalities;
  const hasVariableBounds = variableBounds !== "none";
  let outputText = "";
  outputText += "\\[ \\begin{align*}\n";
  const nonnegativeSlacks = getNonnegativeSlacks(inequalities, equalities, variableBounds);
  outputText += getObjective(logBarriers ? nonnegativeSlacks : []);

  if (hasInequalities || hasEqualities || hasVariableBounds) {
    outputText += "\\text{subject to} \\quad";
  }

  if (hasEqualities) {
    outputText += "& " + A_eq + x + " = " + b_eq + "\\\\\n"
  }

  if (hasInequalities) {
    outputText += "& " + A_ineq + x + " - " + s_A + " = 0 \\\\\n";
    if (inequalities === "lower" || inequalities == "both") outputText += "& " + s_A + " - " + s_lA + " = " + l_A + " \\\\\n";
    if (inequalities === "upper" || inequalities == "both") outputText += "& " + s_A + " + " + s_uA + " = " + u_A + " \\\\\n";
  }

  if (variableBounds !== "none") {
    if (variableBounds === "lower" || variableBounds == "both") outputText += "& " + x + " - " + s_lx + " = " + l_x + " \\\\\n";
    if (variableBounds === "upper" || variableBounds == "both") outputText += "& " + x + " + " + s_ux + " = " + u_x + " \\\\\n";
  }

  if (!logBarriers && nonnegativeSlacks) {
    outputText += "& " + nonnegativeSlacks.join(", ") + " \\geq 0 \\\\\n"
  }

  outputText += "\\end{align*} \\]\n"
  return outputText;
}

function updateProblem() {
  let inequalities = document.querySelector('input[name="inequalities"]:checked').value;
  let handlingInequalities = "Slacks"; // Only one option
  let equalities = document.getElementById("equalities").checked;
  let handlingEqualities = document.querySelector('input[name="handling_equalities"]:checked').value;
  let variableBounds = document.querySelector('input[name="variable_bounds"]:checked').value;

  const hasInequalities = inequalities !== "none";
  const hasEqualities = equalities;
  const hasVariableBounds = variableBounds !== "none";

  let outputText = ""; //<h3>Quadratic Programming Problem</h3>";
  outputText += "<p><strong>Original optimization problem:</strong></p>";
  outputText += getOriginalProblem(inequalities, equalities, variableBounds);

  if (hasInequalities || hasVariableBounds) {
    outputText += "<p><strong>Slacked optimization problem:</strong></p>";
    outputText += getSlackProblem(inequalities, equalities, variableBounds, false);

    outputText += "<p><strong>Slacked optimization problem with log-barriers:</strong></p>";
    outputText += getSlackProblem(inequalities, equalities, variableBounds, true);
  }


  if (wasmModule) {
    const boundsMap = {
      both: wasmModule.Bounds.Both,
      lower: wasmModule.Bounds.Lower,
      upper: wasmModule.Bounds.Upper,
      none: wasmModule.Bounds.None,
    };
    try {
      const settings = new wasmModule.Settings();
      settings.inequalities = boundsMap[inequalities] ?? wasmModule.Bounds.None;
      settings.variableBounds = boundsMap[variableBounds] ?? wasmModule.Bounds.None;
      const lagrangian = wasmModule.getLagrangian(settings);
      outputText += "<p><strong>Lagrangian function:</strong></p>";
      outputText += "\\[ \\begin{align*}" + lagrangian + "\\end{align*} \\]";

      const firstOrder = wasmModule.getFirstOrderOptimalityConditions(settings);
      outputText += "<p><strong>First-order optimality conditions:</strong></p>";
      outputText += "\\[ \\begin{align*}" + firstOrder + "\\end{align*} \\]";

    } catch (error) {
        console.error("Error calling Lagrangian function:", error);
        outputText += "<p>Error generating Lagrangian: " + error.message + "</p>";
    }
  } else {
    outputText += "<p>WASM module not yet initialized. Lagrangian will appear here when ready.</p>";
  }

  outputText += "<h3>TODO list</h3>";
  outputText += "<p><s>1. Optimization problem with slacks</s></p>";
  outputText += "<p><s>2. Optimization problem with barriers</s></p>";
  outputText += "<p><s>3. Lagrangian function</s></p>";
  outputText += "<p>4. First-order optimality conditions</p>";
  outputText += "<p>5. Newton system</p>";
  outputText += "<p>6. Reduction of rows for log-barriers</p>";
  outputText += "<p>7. Reduction of rows for Lagrange multipliers</p>";
  outputText += "<p>8. Augmented system and LU/LDLT solution methods</p>";
  outputText += "<p>9. Reduction to normal equations if possible</p>";
  outputText += "<p>10. Support for equality conditions</p>";
  outputText += "<p>11. Support for direct slacks for inequalities</p>";

  document.getElementById("output").innerHTML = outputText;
  MathJax.typesetPromise().catch((err) => {
    console.error("MathJax error:", err);
  });
}