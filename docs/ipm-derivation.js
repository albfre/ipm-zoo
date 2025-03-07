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
const A_eq = "C"; // Equality constraint matrix
const b_eq = "b"; // Equality constraint right-hand side
const l_x = "l_x"; // Primary variable lower bound
const u_x = "u_x"; // Primary variable upper bound
const s_A = "s"; // Slack variable for inequality constraint
const s_lA = "g"; // Slack variable for inequality constraint lower bound
const s_uA = "t"; // Slack variable for inequality constraint upper bound
const s_lx = "y"; // Slack variable for x lower bound
const s_ux = "z"; // Slack variable for x upper bound
const getObjective = (equalityPenalties, logBarrierVariables = []) => {
  let output = "\\underset{" + x + "}{\\text{minimize}} \\quad & 0.5 " + x + "^T Q " + x + " + c^T" + x;
  if (equalityPenalties) {
    output += "+ 0.5 \\mu^{-1}(" + A_eq + x + " - " + b_eq + ")^T(" + A_eq + x + " - " + b_eq + ")";
  }
  if (logBarrierVariables.length > 0)
  {
    output +=  "\\\\\n & ";
    for (const v of logBarrierVariables) {
      output += "- \\mu e^T \\log(" + v + ")";
    }
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
  outputText += getObjective(false);

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
  console.log("variableBounds")
  console.log(variableBounds)

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
  if (inequalities === "lower" || inequalities === "both") {
    nonnegativeSlacks.push(s_lA);
  }
  if (inequalities === "upper" || inequalities === "both") {
    nonnegativeSlacks.push(s_uA);
  }
  if (variableBounds === "lower" || variableBounds === "both") {
    nonnegativeSlacks.push(s_lx);
  }
  if (variableBounds === "upper" || variableBounds === "both") {
    nonnegativeSlacks.push(s_ux);
  }
  return nonnegativeSlacks;
}

function getSlackProblem(inequalities, equalities, variableBounds, inequalityHandling, equalityHandling, logBarriers) {
  const hasInequalities = inequalities !== "none";
  const hasEqualities = equalities;
  const hasVariableBounds = variableBounds !== "none";
  let outputText = "";
  outputText += "\\[ \\begin{align*}\n";
  const nonnegativeSlacks = getNonnegativeSlacks(inequalities, equalities, variableBounds);
  outputText += getObjective(equalities && equalityHandling === "penalty", logBarriers ? nonnegativeSlacks : []);

  if (hasInequalities || hasEqualities || hasVariableBounds) {
    outputText += "\\text{subject to} \\quad";
  }

  if (hasEqualities && equalityHandling !== "penalty") {
    outputText += "& " + A_eq + x + " = " + b_eq + "\\\\\n"
  }

  if (hasInequalities) {
    const addLower = inequalities === "lower" || inequalities === "both"
    const addUpper = inequalities === "upepr" || inequalities === "both"
    if (inequalityHandling === "slacks") {
      outputText += "& " + A_ineq + x + " - " + s_A + " = 0 \\\\\n";
      if (addLower) outputText += "& " + s_A + " - " + s_lA + " = " + l_A + " \\\\\n";
      if (addUpper) outputText += "& " + s_A + " + " + s_uA + " = " + u_A + " \\\\\n";
    }
    else {
      if (addLower) outputText += "& " + A_ineq + x + " - " + s_lA + " = " + l_A + " \\\\\n";
      if (addUpper) outputText += "& " + A_ineq + x + " + " + s_uA + " = " + u_A + " \\\\\n";
    }
  }

  if (variableBounds !== "none") {
    if (variableBounds === "lower" || variableBounds === "both") outputText += "& " + x + " - " + s_lx + " = " + l_x + " \\\\\n";
    if (variableBounds === "upper" || variableBounds === "both") outputText += "& " + x + " + " + s_ux + " = " + u_x + " \\\\\n";
  }

  if (!logBarriers && nonnegativeSlacks) {
    outputText += "& " + nonnegativeSlacks.join(", ") + " \\geq 0 \\\\\n"
  }

  outputText += "\\end{align*} \\]\n"
  return outputText;
}

function dimZeros(str) {
  const useDimmedZeros = document.getElementById("dim_zeros").checked;
  return useDimmedZeros ? str.replace(/(\D|^)0(\D|$)/g, '$1{\\color{lightgray}0}$2') : str;
}

function countAmpersandsBeforeNewlines(str) {
  const rows = str.split('\\\\');
  rows.pop();
  
  return rows.map(row => {
    return (row.match(/&/g) || []).length;
  });
}

function updateProblem() {
  const inequalities = document.querySelector('input[name="inequalities"]:checked').value;
  const inequalityHandling = document.querySelector('input[name="handling_inequalities"]:checked').value;
  const equalities = document.getElementById("equalities").checked;
  const equalityHandling = document.querySelector('input[name="handling_equalities"]:checked').value;
  const variableBounds = document.querySelector('input[name="variable_bounds"]:checked').value;

  const hasInequalities = inequalities !== "none";
  const hasEqualities = equalities;
  const hasVariableBounds = variableBounds !== "none";

  let outputText = "";
  outputText += "<p><strong>Original optimization problem:</strong></p>";
  outputText += getOriginalProblem(inequalities, equalities, variableBounds);

  if (hasInequalities || hasVariableBounds) {
    outputText += "<p><strong>Slacked optimization problem:</strong></p>";
    outputText += getSlackProblem(inequalities, equalities, variableBounds, inequalityHandling, equalityHandling, false);

    outputText += "<p><strong>Slacked optimization problem with log-barriers:</strong></p>";
    outputText += getSlackProblem(inequalities, equalities, variableBounds, inequalityHandling, equalityHandling, true);
  }


  if (wasmModule) {
    const boundsMap = {
      both: wasmModule.Bounds.Both,
      lower: wasmModule.Bounds.Lower,
      upper: wasmModule.Bounds.Upper,
      none: wasmModule.Bounds.None,
    };
    const equalityHandlingMap ={
      indefinite : wasmModule.EqualityHandling.IndefiniteFactorization,
      slacks : wasmModule.EqualityHandling.Slack,
      penalty : wasmModule.EqualityHandling.PenaltyFunction,
    };
    try {
      const settings = new wasmModule.Settings();
      settings.equalities = equalities;
      settings.inequalities = boundsMap[inequalities] ?? wasmModule.Bounds.None;
      settings.variableBounds = boundsMap[variableBounds] ?? wasmModule.Bounds.None;
      settings.inequalityHandling = inequalityHandling === "slacks" ? wasmModule.InequalityHandling.Slacks : wasmModule.InequalityHandling.SimpleSlacks;
      settings.equalityHandling = equalityHandlingMap[equalityHandling];
      const lagrangian = wasmModule.getLagrangian(settings);
      outputText += "<p><strong>Lagrangian function:</strong></p>";
      outputText += "\\[ \\begin{align*} L =" + lagrangian + "\\end{align*} \\]";

      const firstOrder = wasmModule.getFirstOrderOptimalityConditions(settings);
      outputText += "<p><strong>First-order optimality conditions:</strong></p>";
      outputText += "\\[ \\begin{align*}" + firstOrder + "\\end{align*} \\]";

      const newtonSystem = wasmModule.getNewtonSystem(settings);
      outputText += "<p><strong>Newton system:</strong></p>";
      outputText += "\\[ \\begin{align*} \\nabla^2 L p = -\\nabla L \\end{align*}, \\]";
      outputText += "where"
      let cs = "c".repeat(countAmpersandsBeforeNewlines(newtonSystem.lhs) + 1);
      outputText += "\\[ \\nabla^2 L p = \\left( \\begin{array}{" + cs + "} " + dimZeros(newtonSystem.lhs) + "\\end{array} \\right) "
      outputText += "\\left( \\begin{array}{c} " + newtonSystem.variables + "\\end{array} \\right) \\]";
      outputText += "\\[ = -\\nabla L = \\left( \\begin{array}{c} " + newtonSystem.rhs + "\\end{array} \\right)"
      outputText += "=: \\left( \\begin{array}{c} " + newtonSystem.rhsShorthand + "\\end{array} \\right) \\]";

      const augmentedSystem = wasmModule.getAugmentedSystem(settings);
      outputText += "<p><strong>Augmented system:</strong></p>";
      cs = "c".repeat(countAmpersandsBeforeNewlines(augmentedSystem.lhs) + 1);
      outputText += "\\[ \\left( \\begin{array}{" + cs + "} " + dimZeros(augmentedSystem.lhs) + "\\end{array} \\right) "
      outputText += "\\left( \\begin{array}{c} " + augmentedSystem.variables + "\\end{array} \\right) \\]";
      outputText += "\\[ = \\left( \\begin{array}{ll} " + augmentedSystem.rhs + "\\end{array} \\right) \\]"

      const normalEquation = wasmModule.getNormalEquation(settings);
      outputText += "<p><strong>Normal equation:</strong></p>";
      cs = "c".repeat(countAmpersandsBeforeNewlines(normalEquation.lhs) + 1);
      if (equalities && equalityHandling === "indefinite") {
        outputText += "\\[ \\left( \\begin{array}{" + cs + "} " + dimZeros(normalEquation.lhs) + "\\end{array} \\right) "
        outputText += "\\left( \\begin{array}{c} " + normalEquation.variables + "\\end{array} \\right) \\]";
        outputText += "\\[ = \\left( \\begin{array}{ll} " + normalEquation.rhs + "\\end{array} \\right) \\]"
      }
      else {
        outputText += "\\[ \\left( \\begin{array}{" + cs + "} " + dimZeros(normalEquation.lhs) + "\\end{array} \\right) "
        outputText += normalEquation.variables + " \\]";
        outputText += "\\[ \\begin{array}{ll} =" + normalEquation.rhs + " \\end{array} \\]"
      }


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
  outputText += "<p><s>4. First-order optimality conditions</s></p>";
  outputText += "<p><s>5. Newton system</s></p>";
  outputText += "<p><s>6. Reduction of rows for log-barriers</s></p>";
  outputText += "<p><s>7. Reduction of rows for Lagrange multipliers</s></p>";
  outputText += "<p>8. Expressions for search direction variables in reduced system</p>";
  outputText += "<p><s>9. Augmented system</s></p>";
  outputText += "<p>10. LU/LDLT solution methods</p>";
  outputText += "<p><s>11. Reduction to normal equations if possible</s></p>";
  outputText += "<p><s>12. Support for equalities</s></p>";
  outputText += "<p><s>13. Support for direct slacks for inequalities</s></p>";
  outputText += "<p>14. Support for different equality handling</p>";
  outputText += "<p>15. Correct vector differentiation</p>";

  document.getElementById("output").innerHTML = outputText;
  MathJax.typesetPromise().catch((err) => {
    console.error("MathJax error:", err);
  });
}