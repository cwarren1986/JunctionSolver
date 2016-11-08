using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace TestJunctionSolver
{
    [TestClass]
    public class Tests
    {
        [TestMethod]
        public void Basic()
        {
            int N = 500;
            int M = 500;
            double q = 1.602e-19;
            double eps = 10 * 8.854e-12;
            double thickness = 4.25e-6;
            double Eg = 1.5 * q;
            double Ef0 = -300e-3 * q;
            double ssd = 1e20;
            double tepf = 1e10;
            double V = 1.0;

            var dpl = new System.Collections.Generic.List<JunctionSolver.DefectParameters>();

            dpl.Add(new JunctionSolver.DefectParameters(Ef0, 30e-3 * q, 2 * ssd * q, "Shallow Dopant", "Constant"));

            var deviceParameters = new JunctionSolver.DeviceParameters(N, M, eps, thickness, Eg, Ef0, ssd, tepf, V, 0, dpl);

            JunctionSolver.Device device = new JunctionSolver.Device(deviceParameters);

            JunctionSolver.Utils.Bracket(device);
            
            JunctionSolver.Utils.Solve(device);

            device.CopyResultsToDCArrays();

            JunctionSolver.Utils.Bracket(device, 200, 1e3, 0.030);

            JunctionSolver.Utils.Solve(device, 200, 1e3, 0.030);

            Assert.IsTrue(true);
        }

        [TestMethod]
        public void DriveLevel()
        {
            double Ndl, C0, C1;

            JunctionSolver.DriveLevelCapacitanceProfiling.CalculateDriveLevelDensity(
                10 * JunctionSolver.Constants.VacuumPermittivity,
                new double[] { 30e-3, 40e-3, 50e-3, 60e-3, 70e-3, 80e-3, 90e-3 },
                new double[] { 12.35e-5, 12.30e-5, 12.25e-5, 12.20e-5, 12.15e-5, 12.10e-5, 12.05e-5 },
                out Ndl, out C0, out C1
                );
        }
    }
}
