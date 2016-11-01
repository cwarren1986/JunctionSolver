using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    public static class Utils
    {
        public static double NoumerovStep(double curP, double prevP, double nextQ, double curQ, double prevQ, double stepSizeSquared)
        {
            double nextP = 2.0 * curP - prevP + stepSizeSquared * curQ + (stepSizeSquared / 12) * (nextQ + prevQ - 2 * curQ);
            return nextP;
        }

        public static void NoumerovSolve(Device device)
        {
            for (int i = 0; i < device.NumberOfPositionPoints; i++)
            {
                device.FermiLevel[i] = device.NeutralRegionFermiLevel;
            }
            
            for (int i = 1; i < device.NumberOfPositionPoints - 1; i++)
            {
                device.Potential[i + 1] = 
                    2 * device.Potential[i] 
                    - device.Potential[i - 1]
                    + device.PositionSpacingSquared * device.ElectricFluxConstant * device.ChargeDensity[i];

                /*
                if (device.Potential[i + 1] - device.FermiLevel[i + 1] > device.HalfBandGap)
                {
                    device.FermiLevel[i + 1] = device.Potential[i + 1] - device.HalfBandGap;
                }
                */

                device.ChargeDensity[i + 1] = CalcRho(device.FermiLevel[i + 1] - device.Potential[i + 1], i + 1, device);
                
                device.Potential[i + 1] = NoumerovStep(
                    device.Potential[i],
                    device.Potential[i - 1],
                    device.ElectricFluxConstant * device.ChargeDensity[i + 1],
                    device.ElectricFluxConstant * device.ChargeDensity[i],
                    device.ElectricFluxConstant * device.ChargeDensity[i - 1],
                    device.PositionSpacingSquared
                    );
            }
        }

        public static void NoumerovSolve(Device device, double temperature, double frequency)
        {
            double demarcationEnergy = CalcDemarcationEnergy(device, temperature, frequency);

            /*
            for (int i = 0; i < device.NumberOfPositionPoints; i++)
            {
                device.FermiLevel[i] = device.NeutralRegionFermiLevel;
            }
            */

            for (int i = 1; i < device.NumberOfPositionPoints - 1; i++)
            {
                device.Potential[i + 1] =
                    2 * device.Potential[i]
                    - device.Potential[i - 1]
                    + device.PositionSpacingSquared * device.ElectricFluxConstant * device.ChargeDensity[i];

                /*
                if (device.Potential[i + 1] - device.FermiLevel[i + 1] > device.HalfBandGap)
                {
                    device.FermiLevel[i + 1] = device.Potential[i + 1] - device.HalfBandGap;
                }
                if (device.Potential[i+1] < device.FermiLevel[i+1] + demarcationEnergy)
                {
                    device.ChargeDensity[i + 1] = CalcRho(device.FermiLevel[i + 1] - device.Potential[i + 1], i + 1, device);
                }
                else if (device.DCPotential[i+1] < device.FermiLevel[i+1] + demarcationEnergy)
                {
                    device.ChargeDensity[i + 1] = CalcRho(-demarcationEnergy, i + 1, device);
                }
                else
                {
                    device.ChargeDensity[i + 1] = device.DCChargeDensity[i + 1];
                }
                */
                
                if (device.Potential[i + 1] < device.DCFermiLevel[i + 1] + demarcationEnergy)
                {
                    device.ChargeDensity[i + 1] = CalcRho(device.DCFermiLevel[i + 1] - device.Potential[i + 1], i + 1, device);
                }
                else if (device.DCPotential[i + 1] < device.DCFermiLevel[i + 1] + demarcationEnergy)
                {
                    device.ChargeDensity[i + 1] = CalcRho(-demarcationEnergy, i + 1, device);
                }
                else
                {
                    device.ChargeDensity[i + 1] = device.DCChargeDensity[i + 1];
                }

                device.Potential[i + 1] = NoumerovStep(
                    device.Potential[i],
                    device.Potential[i - 1],
                    device.ElectricFluxConstant * device.ChargeDensity[i + 1],
                    device.ElectricFluxConstant * device.ChargeDensity[i],
                    device.ElectricFluxConstant * device.ChargeDensity[i - 1],
                    device.PositionSpacingSquared
                    );
            }
        }

        public static void RhoTable(Device device)
        {
            for (int j = 0; j < device.NumberOfPositionPoints; j++)
            {
                for (int i = 1; i < device.NumberOfEnergyPoints; i++)
                {
                    device.ChargeDensityTable[j][i] = device.ChargeDensityTable[j][i - 1]
                        + Interp(
                            0.5 * (device.Energy[i] + device.Energy[i - 1]),
                            device.Energy,
                            device.DensityOfStates[j],
                            device.EnergySpacing
                            ) * -device.EnergySpacing;
                }
            }
        }

        public static double CalcRho(double phi, int i, Device device)
        {
            return Interp(phi, device.Energy, device.ChargeDensityTable[i], device.EnergySpacing);
        }

        public static void Calcg(Device device)
        {
            device.DensityOfStates = new double[device.NumberOfPositionPoints][];
            for (int i = 0; i < device.NumberOfPositionPoints; i++)
            {
                device.DensityOfStates[i] = new double[device.NumberOfEnergyPoints];
            }

                foreach (Defect defect in device.DefectList)
            {
                for (int i = 0; i < device.NumberOfPositionPoints; i++)
                {
                    for (int j = 0; j < device.NumberOfEnergyPoints; j++)
                    {
                        device.DensityOfStates[i][j] += defect.DensityOfStates[i][j];
                    }
                }
            }
        }

        public static double[] Gaussian(double[] x, double mu, double sigma)
        {
            double[] result = new double[x.Length];

            for (int i=0; i<x.Length; i++)
            {
                result[i] = Math.Exp(-(x[i] - mu) * (x[i] - mu) / (2 * sigma * sigma)) / (Math.Sqrt(2 * Math.PI) * sigma);
            }

            return result;
        }

        public static void Bracket(Device device)
        {
            device.Message += "Bracketing DC solution...\n";

            device.UpperLimit = device.LowerLimit = 0;

            InitializeDevice(device, 1e-25);
            
            NoumerovSolve(device);

            int numIterations = 0;
            int maxIterations = 20;
            
            if (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage)
            {
                while (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage)
                {
                    device.LowerLimit = device.Potential[0];
                    InitializeDevice(device, 2 * device.LowerLimit);

                    NoumerovSolve(device);

                    numIterations++;
                    
                    if (numIterations >= maxIterations && device.NumRecursion < device.MaxRecursion)
                    {
                        device.NumRecursion++;
                        device.Message += "Having difficulty bracketing.\nThickness is too large.\nReducing thickness by 20%.\n";
                        device.ChangeThickness(device.Thickness * 0.8);
                        Bracket(device);
                        break;
                    }
                }
                device.UpperLimit = 2 * device.Potential[0];
            }
            else
            {
                while (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge > device.Voltage)
                {
                    device.UpperLimit = device.Potential[0];
                    InitializeDevice(device, 0.5 * device.UpperLimit);

                    NoumerovSolve(device);

                    numIterations++;
                    
                    if (numIterations >= maxIterations && device.NumRecursion < device.MaxRecursion)
                    {
                        device.NumRecursion++;
                        device.Message += "Having difficulty bracketing.\nThickness is too large.\nReducing thickness by 20%.\n";
                        device.ChangeThickness(device.Thickness * 0.8);
                        Bracket(device);
                        break;
                    }
                }
                device.LowerLimit = device.Potential[0];
            }
        }

        public static void Bracket(Device device, double temperature, double frequency, double acVoltage)
        {
            device.UpperLimit = device.LowerLimit = 0;

            InitializeDevice(device, device.DCPotential[0]);

            NoumerovSolve(device, temperature, frequency);

            if (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage - acVoltage)
            {
                while (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage - acVoltage)
                {
                    device.LowerLimit = device.Potential[0];
                    InitializeDevice(device, 2 * device.LowerLimit);

                    NoumerovSolve(device, temperature, frequency);
                }
                device.UpperLimit = 2 * device.Potential[0];
            }
            else
            {
                while (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge > device.Voltage - acVoltage)
                {
                    device.UpperLimit = device.Potential[0];
                    InitializeDevice(device, 0.5 * device.UpperLimit);

                    NoumerovSolve(device, temperature, frequency);
                }
                device.LowerLimit = device.Potential[0];
            }
        }

        public static double Calcx0(Device device)
        {
            return Math.Sqrt(device.DielectricConstant / (Constants.ElementaryCharge * 
                Interp(device.NeutralRegionFermiLevel, device.Energy, device.DensityOfStates[0], device.EnergySpacing)));
        }

        public static void Solve(Device device)
        {
            int numIterations = 1;
            double tolerance = 1e-6;
            device.Message += "Finding DC solution...\n";
            device.Message += "     Voltage    Target     Error      Tolerance\n";
            while (Math.Abs(device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge - device.Voltage) > tolerance)
            {
                InitializeDevice(device, 0.5 * (device.UpperLimit + device.LowerLimit));

                NoumerovSolve(device);

                device.Message += String.Format("{0,0:000}  {1,8:f7}  {2,8:f7}  {3,8:e2}  {4,8:e2}\n", 
                    numIterations,
                    device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge, 
                    device.Voltage,
                    Math.Abs(device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge - device.Voltage), 
                    tolerance);

                if (device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge > device.Voltage)
                {
                    device.UpperLimit = device.Potential[0];
                }
                else
                {
                    device.LowerLimit = device.Potential[0];
                }
                if (numIterations >= 50 && device.NumRecursion < device.MaxRecursion)
                {
                    device.NumRecursion++;
                    device.Message += "Maximum number of iterations reached!\nThickness is too large.\nReducing thickness by 20%.\n";
                    device.ChangeThickness(device.Thickness * 0.8);
                    Bracket(device);
                    Solve(device);
                    break;
                }
                numIterations++;
            }

            if (device.ChargeDensity[0] / device.ChargeDensity[device.NumberOfPositionPoints - 1] > 1e-4)
            {
                device.NumRecursion++;
                device.Message += "Thickness is too small.\nIncreasing thickness by 10%.\n";
                device.ChangeThickness(device.Thickness * 1.1);
                Bracket(device);
                Solve(device);
            }
        }

        public static void Solve(Device device, double temperature, double frequency, double acVoltage)
        {
            double tolerance = 1e-6;

            while (Math.Abs(device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge - device.Voltage + acVoltage) > tolerance)
            {
                InitializeDevice(device, 0.5 * (device.UpperLimit + device.LowerLimit));

                NoumerovSolve(device, temperature, frequency);
                
                if (device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge > device.Voltage - acVoltage)
                {
                    device.UpperLimit = device.Potential[0];
                }
                else
                {
                    device.LowerLimit = device.Potential[0];
                }
            }
        }

        public static double Interp(double x, double[] xarray, double[] yarray, double dx)
        {
            double result = 0;
            int i = (int)Math.Floor((x - xarray[0])/ dx);
            //if (i <= 0) { return yarray[0]; }
            if (i >= yarray.Length - 1) { return yarray[yarray.Length - 1]; }
            result = yarray[i] + (x - xarray[i]) * (yarray[i + 1] - yarray[i]) / (xarray[i + 1] - xarray[i]);
            return result;
        }

        public static void InitializeDevice(Device device, double initialPotential)
        {
            device.Potential[0] = initialPotential;
            device.Potential[1] = device.Potential[0] * Math.Exp(device.PositionSpacing / device.X0);
            device.ChargeDensity[0] = CalcRho(device.FermiLevel[0] - device.Potential[0], 0, device);
            device.ChargeDensity[1] = CalcRho(device.FermiLevel[1] - device.Potential[1], 1, device);
        }

        public static double CalcDemarcationEnergy(Device device, double temperature, double frequency)
        {
            return Constants.BoltzmannConstant * temperature * Math.Log(device.ThermalEmissionPrefactor * temperature * temperature / frequency);
        }

        public static double CalcCapacitance(Device device)
        {
            double dQ = 0;
            double dV = (device.Potential[device.NumberOfPositionPoints - 1] - device.DCPotential[device.NumberOfPositionPoints - 1]) / Constants.ElementaryCharge;
            for (int i = 0; i < device.NumberOfPositionPoints; i++)
            {
                dQ += device.ChargeDensity[i] - device.DCChargeDensity[i];
            }
            return dQ / dV * device.PositionSpacing;
        }

        public static void CalcDriveLevelDensity(Device device, double[] driveLevelValues, double[] capacitanceValues, 
            out double Ndl, out double C0, out double C1)
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
            Ndl = -C0 * C0 * C0 / (2 * Constants.ElementaryCharge * device.DielectricConstant * C1);
        }

    }
}
