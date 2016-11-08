using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    /// <summary>
    /// Tools for doing drive level capacitance profiling.
    /// </summary>
    public static class DriveLevelCapacitanceProfiling
    {

        #region [Public Static Methods]

        /// <summary>
        /// Calculates the drive level density using a Levenberg-Marquardt non-linear least squares solver.
        /// </summary>
        /// <param name="dielectricConstant">The dielectric constant (in F/m) of the device.</param>
        /// <param name="acVoltageValues">An array of AC voltages (in V) applied to the device.</param>
        /// <param name="capacitanceValues">An array of capacitance values (in F/m^2) corresponding to the AC voltages.</param>
        /// <param name="Ndl">The drive level density in (/m^3).</param>
        /// <param name="C0">The 1st coefficient in Taylor expansion of the capacitance (in F/m^2).</param>
        /// <param name="C1">The 2nd coefficient in Taylor expansion of the capacitance (in F/m^2/V).</param>
        public static void CalculateDriveLevelDensity(
            double dielectricConstant, 
            double[] acVoltageValues, 
            double[] capacitanceValues,
            out double Ndl, 
            out double C0, 
            out double C1
            )
        {
            // Create an instance of the Levenberg-Marquardt non-linear least squares solver.
            LMDotNet.LMSolver lmsolver = new LMDotNet.LMSolver();

            // Fit the data using the 3rd order Taylor expansion model.
            var fit = lmsolver.FitCurve(
                (V, C) => C[0] - C[1] * V + 2 * C[1] * C[1] / C[0] * V * V - 5 * C[1] * C[1] * C[1] / (C[0] * C[0]) * V * V * V,
                new[] { 1e-5, -1e-5 }, // Guess that C0 = 1e-5 and C1 = -1e-5.
                acVoltageValues,
                capacitanceValues
                );

            // Extract C0 and C1.
            C0 = fit.OptimizedParameters[0];
            C1 = fit.OptimizedParameters[1];

            // Calculate the drive level density.
            Ndl = -C0 * C0 * C0 / (2 * Constants.ElementaryCharge * dielectricConstant * C1);
        }

        #endregion [Public Static Methods]
    }
}
