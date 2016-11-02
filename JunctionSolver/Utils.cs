using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    /// <summary>
    /// Utilities for solving Poisson's equation.
    /// </summary>
    internal static class Utils
    {
        /// <summary>
        /// Calculate the next value of P on the grid using Numerov method, where P'' = Q.
        /// </summary>
        /// <param name="curP">P at position i</param>
        /// <param name="prevP">P at position i-1</param>
        /// <param name="nextQ">Q at position i+1</param>
        /// <param name="curQ">Q at position i</param>
        /// <param name="prevQ">Q at position i-1</param>
        /// <param name="stepSizeSquared">The squared grid spacing.</param>
        /// <returns>P at position i+1</returns>
        private static double NumerovStep(double curP, double prevP, double nextQ, double curQ, double prevQ, double stepSizeSquared)
        {
            double nextP = 2.0 * curP - prevP + stepSizeSquared * curQ + (stepSizeSquared / 12) * (nextQ + prevQ - 2 * curQ);
            return nextP;
        }

        /// <summary>
        /// Calculate the next potential value using a Numerov step on Poisson's equation.
        /// </summary>
        /// <param name="device">The device being solved.</param>
        /// <param name="currentIndex">The current index of the position grid.</param>
        /// <returns></returns>
        private static double PoissonNumerovStep(Device device, int currentIndex)
        {
            // The charge density needs to be multiplied by the electric flux constant (q/ε)
            // to get the equation into the form required by the general Numerov method (P'' = Q)
            return NumerovStep(
                    device.Potential[currentIndex],
                    device.Potential[currentIndex - 1],
                    device.ElectricFluxConstant * device.ChargeDensity[currentIndex + 1],
                    device.ElectricFluxConstant * device.ChargeDensity[currentIndex],
                    device.ElectricFluxConstant * device.ChargeDensity[currentIndex - 1],
                    device.PositionSpacingSquared
                    );
        }

        /// <summary>
        /// Estimates the potential at currentIndex + 1.
        /// </summary>
        /// <remarks>
        /// Importantly, this does not use any information about the charge density at currentIndex + 1.
        /// Thus, we can use this to estimate the next potential value, which will then allow us to 
        /// estimate the charge density at currentIndex + 1, which in turn allows us to calculate a more
        /// accurate value of the potential at currentIndex + 1.
        /// </remarks>
        /// <param name="device">The device object being solved.</param>
        /// <param name="currentIndex">The currentIndex on the position grid.</param>
        /// <returns>An estimate of the potential at currentIndex + 1.</returns>
        private static double EstimateNextPotential(Device device, int currentIndex)
        {
            return 2 * device.Potential[currentIndex]
                    - device.Potential[currentIndex - 1]
                    + device.PositionSpacingSquared * device.ElectricFluxConstant * device.ChargeDensity[currentIndex];
        }

        /// <summary>
        /// Solves Poisson's equation using the Numerov method.
        /// </summary>
        /// <remarks>
        /// Only makes a single pass over the grid. Solution is entirely determine by the initial values 
        /// of the potential (at the back of the device) and the density of states. Thus, to statisfy the 
        /// boundary conditions at the interface, this method will need to be part of a larger scheme, e.g.,
        /// the shooting method.
        /// </remarks>
        /// <param name="device">The device object to solve.</param>
        private static void NumerovSolve(Device device)
        {
            /* Force the fermi level to be flat. Not needed right now because the Fermi level is never changed.
            for (int i = 0; i < device.NumberOfPositionPoints; i++)
            {
                device.FermiLevel[i] = device.NeutralRegionFermiLevel;
            }
            */
            
            // Loop over the position grid starting at the second position (i=1).
            for (int i = 1; i < device.NumberOfPositionPoints - 1; i++)
            {
                // Estimate the next potential value.
                device.Potential[i + 1] = EstimateNextPotential(device, i);

                /* Adjust the Fermi level so it doesn't cross mid-gap. Not currently used.
                if (device.Potential[i + 1] - device.FermiLevel[i + 1] > device.HalfBandGap)
                {
                    device.FermiLevel[i + 1] = device.Potential[i + 1] - device.HalfBandGap;
                }
                */

                // Calculate the next charge density value based on the estimate.
                device.ChargeDensity[i + 1] = CalculateChargeDensity(device.FermiLevel[i + 1] - device.Potential[i + 1], i + 1, device);

                // Recalculate the next potential value using a Numerov step.
                device.Potential[i + 1] = PoissonNumerovStep(device, i);
            }
        }

        private static void NoumerovSolve(Device device, double temperature, double frequency)
        {
            double demarcationEnergy = CalcDemarcationEnergy(device, temperature, frequency);

            /* Force the fermi level to be flat. Not needed right now because the Fermi level is never changed.
            for (int i = 0; i < device.NumberOfPositionPoints; i++)
            {
                device.FermiLevel[i] = device.NeutralRegionFermiLevel;
            }
            */

            for (int i = 1; i < device.NumberOfPositionPoints - 1; i++)
            {
                // Estimate the next AC potential value.
                device.Potential[i + 1] = EstimateNextPotential(device, i);

                /* Adjust the Fermi level so it doesn't cross mid-gap. Not currently used.
                if (device.Potential[i + 1] - device.FermiLevel[i + 1] > device.HalfBandGap)
                {
                    device.FermiLevel[i + 1] = device.Potential[i + 1] - device.HalfBandGap;
                }
                */

                // Calculate next AC charge density accounting for whether or not the charge can respond.
                if (device.Potential[i + 1] < device.DCFermiLevel[i + 1] + demarcationEnergy)
                {
                    // The AC potential is within a demarcation energy of the Fermi level-- everything can respond.
                    // Calculate the next charge density as normal.
                    device.ChargeDensity[i + 1] = CalculateChargeDensity(device.DCFermiLevel[i + 1] - device.Potential[i + 1], i + 1, device);
                }
                else if (device.DCPotential[i + 1] < device.DCFermiLevel[i + 1] + demarcationEnergy)
                {
                    // The AC potenial is not within a demarcation energy of the fermi level (b/c the above condition failed),
                    // but the DC potenital still is. The release of charge is emission limited, therefore the charge density
                    // is simply detemined by the demarcation energy.
                    device.ChargeDensity[i + 1] = CalculateChargeDensity(-demarcationEnergy, i + 1, device);
                }
                else
                {
                    // Neither the AC or DC potential is within a demarcation energy of the Fermi-- nothing can respond.
                    // The charge density is frozen in the DC configuration. 
                    device.ChargeDensity[i + 1] = device.DCChargeDensity[i + 1];
                }

                // Recalculate the next potential value using a Numerov step.
                device.Potential[i + 1] = PoissonNumerovStep(device, i);
            }
        }

        internal static void RhoTable(Device device)
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

        private static double CalculateChargeDensity(double phi, int i, Device device)
        {
            return Interp(phi, device.Energy, device.ChargeDensityTable[i], device.EnergySpacing);
        }

        internal static void CalculateDensityOfStates(Device device)
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

        internal static double[] Gaussian(double[] x, double mu, double sigma)
        {
            double[] result = new double[x.Length];

            for (int i=0; i<x.Length; i++)
            {
                result[i] = Math.Exp(-(x[i] - mu) * (x[i] - mu) / (2 * sigma * sigma)) / (Math.Sqrt(2 * Math.PI) * sigma);
            }

            return result;
        }

        internal static void Bracket(Device device)
        {
            device.Message += "Bracketing DC solution...\n";

            device.UpperLimit = device.LowerLimit = 0;

            InitializeDevice(device, 1e-25);
            
            NumerovSolve(device);

            int numIterations = 0;
            int maxIterations = 20;
            
            if (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage)
            {
                while (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage)
                {
                    device.LowerLimit = device.Potential[0];
                    InitializeDevice(device, 2 * device.LowerLimit);

                    NumerovSolve(device);

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

                    NumerovSolve(device);

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

        internal static void Bracket(Device device, double temperature, double frequency, double acVoltage)
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

        internal static double Calcx0(Device device)
        {
            return Math.Sqrt(device.DielectricConstant / (Constants.ElementaryCharge * 
                Interp(device.NeutralRegionFermiLevel, device.Energy, device.DensityOfStates[0], device.EnergySpacing)));
        }

        internal static void Solve(Device device)
        {
            int numIterations = 1;
            double tolerance = 1e-6;
            device.Message += "Finding DC solution...\n";
            device.Message += "     Voltage    Target     Error      Tolerance\n";
            while (Math.Abs(device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge - device.Voltage) > tolerance)
            {
                InitializeDevice(device, 0.5 * (device.UpperLimit + device.LowerLimit));

                NumerovSolve(device);

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

        internal static void Solve(Device device, double temperature, double frequency, double acVoltage)
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

        private static double Interp(double x, double[] xarray, double[] yarray, double dx)
        {
            double result = 0;
            int i = (int)Math.Floor((x - xarray[0])/ dx);
            //if (i <= 0) { return yarray[0]; }
            if (i >= yarray.Length - 1) { return yarray[yarray.Length - 1]; }
            result = yarray[i] + (x - xarray[i]) * (yarray[i + 1] - yarray[i]) / (xarray[i + 1] - xarray[i]);
            return result;
        }

        private static void InitializeDevice(Device device, double initialPotential)
        {
            device.Potential[0] = initialPotential;
            device.Potential[1] = device.Potential[0] * Math.Exp(device.PositionSpacing / device.X0);
            device.ChargeDensity[0] = CalculateChargeDensity(device.FermiLevel[0] - device.Potential[0], 0, device);
            device.ChargeDensity[1] = CalculateChargeDensity(device.FermiLevel[1] - device.Potential[1], 1, device);
        }

        internal static double CalcDemarcationEnergy(Device device, double temperature, double frequency)
        {
            return Constants.BoltzmannConstant * temperature * Math.Log(device.ThermalEmissionPrefactor * temperature * temperature / frequency);
        }

        internal static double CalcCapacitance(Device device)
        {
            double dV = (device.Potential[device.NumberOfPositionPoints - 1] 
                - device.DCPotential[device.NumberOfPositionPoints - 1]) 
                / Constants.ElementaryCharge;

            double dQ = 0;
            for (int i = 0; i < device.NumberOfPositionPoints; i++)
            {
                dQ += device.ChargeDensity[i] - device.DCChargeDensity[i];
            }
            return dQ / dV * device.PositionSpacing;
        }
    }
}
