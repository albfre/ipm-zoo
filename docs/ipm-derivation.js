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
const l_A_ineq = "l_{A}"; // Inequality constraint lower bound
const u_A_ineq = "u_{A}"; // Inequality constraint upper bound
const A_eq = "C"; // Equality constraint matrix
const b_eq = "d"; // Equality constraint right-hand side
const p_eq = "p"; // Equality constraint regularization
const delta_eq = "\\delta"; // Equality constraint regularization factor
const l_x = "l_x"; // Primary variable lower bound
const u_x = "u_x"; // Primary variable upper bound
const s_A_ineq = "s"; // Slack variable for inequality constraint
const s_lA = "g"; // Slack variable for inequality constraint lower bound
const s_uA = "h"; // Slack variable for inequality constraint upper bound
const s_lx = "y"; // Slack variable for x lower bound
const s_ux = "z"; // Slack variable for x upper bound
const s_A_eq = "t"; // Slack variable for equality constraints
const s_lC = "v"; // Slack variable for equality constraint lower bound
const s_uC = "w"; // Slack variable for equality constraint upper bound

const getObjective = (equalityPenalties = false, equalityRegularization = false, logBarrierVariables = []) => {
  let output = "\\text{minimize} \\quad & 0.5 " + x + "^T Q " + x + " + c^T" + x;
  if (equalityPenalties) {
    output += "+ 0.5 \\mu^{-1}(" + A_eq + x + " - " + b_eq + ")^T(" + A_eq + x + " - " + b_eq + ")";
  }
  if (equalityRegularization) {
    output += "+ 0.5 " + p_eq + "^T" + p_eq;

  }
  if (logBarrierVariables.length > 0)
  {
    output +=  " \\\\\n & ";
    for (let i = 0; i < logBarrierVariables.length; ++i) {
      if (logBarrierVariables.length > 4 && i === logBarrierVariables.length / 2) {
        output += " \\\\\n & "
      }
      output += "- \\mu e^T \\log(" + logBarrierVariables[i] + ")";
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
  outputText += getObjective();

  if (hasInequalities || hasEqualities || hasVariableBounds) {
    outputText += "\\text{subject to} \\quad";
  }

  if (hasInequalities) {
    if (inequalities === "lower") outputText += "& " + A_ineq + x + " \\geq " + l_A_ineq + " \\\\\n";
    if (inequalities === "upper") outputText += "& " + A_ineq + x + " \\leq " + u_A_ineq + " \\\\\n";
    if (inequalities === "both") outputText += "& " + l_A_ineq + " \\leq " + A_ineq + x + " \\leq " + u_A_ineq + " \\\\\n";
  }

  if (hasEqualities) {
    outputText += "& " + A_eq + x + " = " + b_eq + "\\\\\n"
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
  if (inequalities === "lower" || inequalities === "both") {
    nonnegativeSlacks.push(s_lA);
  }
  if (inequalities === "upper" || inequalities === "both") {
    nonnegativeSlacks.push(s_uA);
  }
  if (equalities) {
    nonnegativeSlacks.push(s_lC);
    nonnegativeSlacks.push(s_uC);
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
  const hasSlackedEqualities = hasEqualities && (equalityHandling === "slacks" || equalityHandling === "simple_slacks");
  const nonnegativeSlacks = getNonnegativeSlacks(inequalities, hasSlackedEqualities, variableBounds);
  const equalityPenalties = equalities && equalityHandling === "penalty";
  const equalityRegularization = equalities && equalityHandling === "regularization";
  outputText += getObjective(equalityPenalties, equalityRegularization, logBarriers ? nonnegativeSlacks : []);

  if (hasInequalities || hasEqualities || hasVariableBounds) {
    outputText += "\\text{subject to} \\quad";
  }

  if (hasInequalities) {
    const addLower = inequalities === "lower" || inequalities === "both"
    const addUpper = inequalities === "upepr" || inequalities === "both"
    if (inequalityHandling === "slacks") {
      outputText += "& " + A_ineq + x + " - " + s_A_ineq + " = 0 \\\\\n";
      if (addLower) outputText += "& " + s_A_ineq + " - " + s_lA + " = " + l_A_ineq + " \\\\\n";
      if (addUpper) outputText += "& " + s_A_ineq + " + " + s_uA + " = " + u_A_ineq + " \\\\\n";
    }
    else {
      if (addLower) outputText += "& " + A_ineq + x + " - " + s_lA + " = " + l_A_ineq + " \\\\\n";
      if (addUpper) outputText += "& " + A_ineq + x + " + " + s_uA + " = " + u_A_ineq + " \\\\\n";
    }
  }

  if (hasEqualities && equalityHandling !== "penalty") {
    if (equalityHandling === "slacks") {
      outputText += "& " + A_eq + x + " - " + s_A_eq + " = 0 \\\\\n";
      outputText += "& " + s_A_eq + " - " + s_lC + " = " + b_eq + " \\\\\n";
      outputText += "& " + s_A_eq + " + " + s_uC + " = " + b_eq + " \\\\\n";
    }
    else if (equalityHandling == "simple_slacks") {
      outputText += "& " + A_eq + x + " - " + s_lC + " = " + b_eq + " \\\\\n";
      outputText += "& " + A_eq + x + " + " + s_uC + " = " + b_eq + " \\\\\n";
    }
    else if (equalityHandling == "regularization") {
      outputText += "& " + A_eq + x + " + " + delta_eq + " " + p_eq + " = " + b_eq + " \\\\\n";
    }
    else {
      outputText += "& " + A_eq + x + " = " + b_eq + " \\\\\n";
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
  return rows.length > 0 ? (rows[0].match(/&/g) || []).length : 0;
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

  const hasSlackedEqualities = hasEqualities && (equalityHandling === "slacks" || equalityHandling === "simple_slacks");
  if (hasInequalities || hasVariableBounds || hasSlackedEqualities) {
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
      none: wasmModule.EqualityHandling.None,
      slacks : wasmModule.EqualityHandling.Slacks,
      simple_slacks: wasmModule.EqualityHandling.SimpleSlacks,
      penalty : wasmModule.EqualityHandling.PenaltyFunction,
      regularization : wasmModule.EqualityHandling.Regularization,
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
      outputText += "\\[ \\begin{align*}\n L =" + lagrangian + "\\end{align*} \\]";

      const firstOrder = wasmModule.getFirstOrderOptimalityConditions(settings);
      outputText += "<p><strong>First-order optimality conditions:</strong></p>";
      outputText += "\\[ \\begin{align*}\n " + firstOrder + "\\end{align*} \\]";

      const newtonSystems = wasmModule.getNewtonSystems(settings);
      const newtonSystem = newtonSystems.full;
      outputText += "<p><strong>Newton system:</strong></p>";
      outputText += "\\[ \\begin{align*}\n \\nabla^2 L p = -\\nabla L \\end{align*}, \\]";
      outputText += "or"
      let cs = "c".repeat(countAmpersandsBeforeNewlines(newtonSystem.lhs) + 1);
      outputText += "\\[ \\left( \\begin{array}{" + cs + "}\n " + dimZeros(newtonSystem.lhs) + "\\end{array} \\right) "
      outputText += "\\left( \\begin{array}{c}\n " + newtonSystem.variables + "\\end{array} \\right) \\]";
      outputText += "\\[ = \\left( \\begin{array}{c}\n " + newtonSystem.rhs + "\\end{array} \\right)"
      outputText += "=: \\left( \\begin{array}{c}\n " + newtonSystem.rhsShorthand + "\\end{array} \\right) \\]";

      const augmentedSystem = newtonSystems.augmented;
      outputText += "<p><strong>Augmented system:</strong></p>";
      cs = "c".repeat(countAmpersandsBeforeNewlines(augmentedSystem.lhs) + 1);
      outputText += "\\[ \\left( \\begin{array}{" + cs + "}\n " + dimZeros(augmentedSystem.lhs) + "\\end{array} \\right) "
      outputText += "\\left( \\begin{array}{c}\n " + augmentedSystem.variables + "\\end{array} \\right) \\]";
      outputText += "\\[ = \\left( \\begin{array}{l}\n " + augmentedSystem.rhs + "\\end{array} \\right) \\]"
      outputText += "where"
      outputText += "\\[ \\begin{align*}\n " + augmentedSystem.deltaDefinitions + "\\end{align*} \\]"

      const normalEquations = newtonSystems.normal;
      outputText += "<p><strong>Normal equations:</strong></p>";
      const numNormalEquationVariables = countAmpersandsBeforeNewlines(normalEquations.lhs) + 1;
      cs = "c".repeat(numNormalEquationVariables);
      outputText += "\\[ \\left( \\begin{array}{" + cs + "}\n " + dimZeros(normalEquations.lhs) + "\\end{array} \\right) "
      if (numNormalEquationVariables > 1) {
        outputText += "\\left( \\begin{array}{c}\n " + normalEquations.variables + "\\end{array} \\right) \\]";
        outputText += "\\[ = \\left( \\begin{array}{l}\n " + normalEquations.rhs + "\\end{array} \\right) \\]"
      }
      else {
        outputText += normalEquations.variables + " \\]";
        outputText += "\\[ \\begin{array}{c}\n =" + normalEquations.rhs + " \\end{array} \\]"
      }
      outputText += "where"
      outputText += "\\[ \\begin{align*}\n " + normalEquations.deltaDefinitions + "\\end{align*} \\]"
    } catch (error) {
      console.error("Error calling Lagrangian function:", error);
      outputText += "<p>Error generating Lagrangian: " + error.message + "</p>";
    }
  } else {
    outputText += "<p>WASM module not yet initialized. Lagrangian will appear here when ready.</p>";
  }

  document.getElementById("output").innerHTML = outputText;
  MathJax.typesetPromise().catch((err) => {
    console.error("MathJax error:", err);
  });
}