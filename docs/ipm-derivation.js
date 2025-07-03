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

  if (wasmModule) {
    const boundsMap = {
      both: wasmModule.Bounds.Both,
      lower: wasmModule.Bounds.Lower,
      upper: wasmModule.Bounds.Upper,
      none: wasmModule.Bounds.None,
    };
    const inequalityHandlingMap = {
      slacks : wasmModule.InequalityHandling.Slacks,
      slacked_slacks : wasmModule.InequalityHandling.SlackedSlacks,
      naive_slacks: wasmModule.InequalityHandling.NaiveSlacks,
    };
    const equalityHandlingMap = {
      none: wasmModule.EqualityHandling.None,
      slacks : wasmModule.EqualityHandling.Slacks,
      slacked_slacks : wasmModule.EqualityHandling.SlackedSlacks,
      naive_slacks: wasmModule.EqualityHandling.NaiveSlacks,
      penalty : wasmModule.EqualityHandling.PenaltyFunction,
      regularization : wasmModule.EqualityHandling.Regularization,
    };
    try {
      const settings = new wasmModule.Settings();
      settings.equalities = equalities;
      settings.inequalities = boundsMap[inequalities] ?? wasmModule.Bounds.None;
      settings.variableBounds = boundsMap[variableBounds] ?? wasmModule.Bounds.None;
      settings.inequalityHandling = inequalityHandlingMap[inequalityHandling];
      settings.equalityHandling = equalityHandlingMap[equalityHandling];

      outputText += "<p><strong>Original optimization problem:</strong></p>";
      const original = wasmModule.getOptimizationProblem(settings, wasmModule.OptimizationProblemType.Original);
      outputText += "\\[ \\begin{align*}\n " + original + "\\end{align*} \\]";

      const hasSlackedEqualities = hasEqualities && (equalityHandling === "slacks" || equalityHandling === "slacked_slacks" || equalityHandling === "naive_slacks");
      if (hasInequalities || hasVariableBounds || hasSlackedEqualities) {
        outputText += "<p><strong>Slacked optimization problem:</strong></p>";
        const slacked = wasmModule.getOptimizationProblem(settings, wasmModule.OptimizationProblemType.Slacked);
        outputText += "\\[ \\begin{align*}\n " + slacked + "\\end{align*} \\]";

        outputText += "<p><strong>Slacked optimization problem with log-barriers:</strong></p>";
        const barrier = wasmModule.getOptimizationProblem(settings, wasmModule.OptimizationProblemType.SlackedWithBarriers);
        outputText += "\\[ \\begin{align*}\n " + barrier + "\\end{align*} \\]";
      }

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
      if (augmentedSystem.variables !== "") {
        outputText += "<p><strong>Augmented system:</strong></p>";
        cs = "c".repeat(countAmpersandsBeforeNewlines(augmentedSystem.lhs) + 1);
        outputText += "\\[ \\left( \\begin{array}{" + cs + "}\n " + dimZeros(augmentedSystem.lhs) + "\\end{array} \\right) "
        outputText += "\\left( \\begin{array}{c}\n " + augmentedSystem.variables + "\\end{array} \\right) \\]";
        outputText += "\\[ = \\left( \\begin{array}{l}\n " + augmentedSystem.rhs + "\\end{array} \\right) \\]"
        outputText += "where"
        outputText += "\\[ \\begin{align*}\n " + augmentedSystem.deltaDefinitions + "\\end{align*} \\]"
      }

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