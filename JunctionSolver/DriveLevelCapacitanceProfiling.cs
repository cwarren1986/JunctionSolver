using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    public static class DriveLevelCapacitanceProfiling
    {
        public static void CalculateDriveLevelDensity(
            double dielectricConstant, 
            double[] driveLevelValues, 
            double[] capacitanceValues,
            out double Ndl, 
            out double C0, 
            out double C1
            )
        {
            LMDotNet.LMSolver lmsolver = new LMDotNet.LMSolver();
            var fit = lmsolver.FitCurve(
                (V, C) => C[0] - C[1] * V + 2 * C[1] * C[1] / C[0] * V * V - 5 * C[1] * C[1] * C[1] / (C[0] * C[0]) * V * V * V,
                new[] { 1e-5, -1e-5 },
                driveLevelValues,
                capacitanceValues
                );
            C0 = fit.OptimizedParameters[0];
            C1 = fit.OptimizedParameters[1];
            Ndl = -C0 * C0 * C0 / (2 * Constants.ElementaryCharge * dielectricConstant * C1);
        }
    }
}
