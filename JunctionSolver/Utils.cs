using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    /// <summary>
    /// Utilities for solving Poisson's equation in a device.
    /// </summary>
    internal static class Utils
    {

        #region [Private Methods]

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

        /// <summary>
        /// Overload of NumerovSolve that calculates the AC solution.
        /// </summary>
        /// <param name="device">The device being solved.</param>
        /// <param name="temperature">The temperature (in K) of the device.</param>
        /// <param name="frequency">The frequency (in Hz) of the AC perturbation.</param>
        private static void NoumerovSolve(Device device, double temperature, double frequency)
        {
            double demarcationEnergy = CalculateDemarcationEnergy(device, temperature, frequency);

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

        /// <summary>
        /// Calculates the charge density.
        /// </summary>
        /// <param name="phi">The potential (in J).</param>
        /// <param name="i">The position index.</param>
        /// <param name="device">The device being solved.</param>
        /// <returns>The charge density in (C/m^3).</returns>
        private static double CalculateChargeDensity(double phi, int i, Device device)
        {
            // The charge density integrals are precalculated and stored in ChargeDensityIntegralTable.
            // To calculate the charge density integral simply interpolate.
            return Interpolate(phi, device.Energy, device.ChargeDensityIntegralTable[i], device.EnergySpacing);
        }

        /// <summary>
        /// A simple linear interpolation routine.
        /// </summary>
        /// <param name="x">The x-value to be interpolated to.</param>
        /// <param name="xarray">The array x-values.</param>
        /// <param name="yarray">The array of y-values</param>
        /// <param name="dx">The spacing of xarray.</param>
        /// <returns>The value of y at x.</returns>
        private static double Interpolate(double x, double[] xarray, double[] yarray, double dx)
        {
            // Calculate the last index i of xarray such that xarray[i] < x
            int i = (int)Math.Floor((x - xarray[0]) / dx);

            // If i is outside the allowed range, set result to the appropriate end of yarray.
            if (i <= 0) { return yarray[0]; }
            if (i >= yarray.Length - 1) { return yarray[yarray.Length - 1]; }

            // If i is legal, interpolate and return.
            return yarray[i] + (x - xarray[i]) * (yarray[i + 1] - yarray[i]) / (xarray[i + 1] - xarray[i]);
        }

        /// <summary>
        /// Initializes the values of the potential and the charge density at the basck of the device.
        /// </summary>
        /// <param name="device">The device to be solved.</param>
        /// <param name="initialPotential">The value to which Potential[0] will be set.</param>
        private static void InitializeBoundary(Device device, double initialPotential)
        {
            // Set Potential[0].
            device.Potential[0] = initialPotential;

            // Calculate the next potential value assuming exponential bands following the constant DOS formula.
            device.Potential[1] = device.Potential[0] * Math.Exp(device.PositionSpacing / device.X0);

            // Calculate the first and second charge density values consistent with the potential.
            device.ChargeDensity[0] = CalculateChargeDensity(device.FermiLevel[0] - device.Potential[0], 0, device);
            device.ChargeDensity[1] = CalculateChargeDensity(device.FermiLevel[1] - device.Potential[1], 1, device);
        }

        #endregion [Private Methods]

        #region [Internal Methods]

        /// <summary>
        /// Precalculates the charge density integrals.
        /// </summary>
        /// <param name="device">The device to be solved.</param>
        internal static void CalculateChargeDensityIntegralTable(Device device)
        {
            // Loop over the position grid.
            for (int j = 0; j < device.NumberOfPositionPoints; j++)
            {
                // Loop over the energy grid.
                for (int i = 1; i < device.NumberOfEnergyPoints; i++)
                {
                    // Calculate the current charge density integral from the previous
                    // one by adding on the next incremental piece of the full integral.
                    device.ChargeDensityIntegralTable[j][i] = device.ChargeDensityIntegralTable[j][i - 1]
                        + Interpolate(
                            0.5 * (device.Energy[i] + device.Energy[i - 1]),
                            device.Energy,
                            device.DensityOfStates[j],
                            device.EnergySpacing
                            ) * -device.EnergySpacing;
                }
            }
        }

        /// <summary>
        /// Calculates the density of states.
        /// </summary>
        /// <param name="device"></param>
        internal static void CalculateDensityOfStates(Device device)
        {
            // Initialize the density of states array.
            device.DensityOfStates = new double[device.NumberOfPositionPoints][];
            for (int i = 0; i < device.NumberOfPositionPoints; i++)
            {
                device.DensityOfStates[i] = new double[device.NumberOfEnergyPoints];
            }

            // Loop over DefectList.
            foreach (Defect defect in device.DefectList)
            {
                // Loop over the position grid.
                for (int i = 0; i < device.NumberOfPositionPoints; i++)
                {
                    // Loop over the energy grid.
                    for (int j = 0; j < device.NumberOfEnergyPoints; j++)
                    {
                        // Add the defect to the density of states of the device. 
                        device.DensityOfStates[i][j] += defect.DensityOfStates[i][j];
                    }
                }
            }
        }
        
        /// <summary>
        /// Calculates an array of values following a gaussian distribution.
        /// </summary>
        /// <param name="x">An array of values at which to evaluate the Gaussian.</param>
        /// <param name="mu">The mean of the Gaussian.</param>
        /// <param name="sigma">The standard deviation of the Gaussian.</param>
        /// <param name="magnitude">The magnitude of the gaussian, i.e., what it integrates to.</param>
        /// <returns>An array of values following a gaussian distribution. Integrates to magnitude.</returns>
        internal static double[] Gaussian(double[] x, double mu, double sigma, double magnitude)
        {
            // Initialize the result array.
            double[] result = new double[x.Length];

            // Loop over x.
            for (int i = 0; i < x.Length; i++)
            {
                // Calculate the Gaussian
                result[i] = magnitude * Math.Exp(-(x[i] - mu) * (x[i] - mu) / (2 * sigma * sigma)) / (Math.Sqrt(2 * Math.PI) * sigma);
            }

            return result;
        }

        /// <summary>
        /// Brackets the initial condition for the DC case.
        /// </summary>
        /// <param name="device">The device to be solved.</param>
        internal static void Bracket(Device device)
        {
            device.Message += "Bracketing DC solution...\n";

            // Set the brackets to zero.
            device.UpperLimit = device.LowerLimit = 0;

            // Initialize the potenial and charge at the back of the device.
            InitializeBoundary(device, 1e-25);
            
            // Solve for the potential and the charge density profiles.
            NumerovSolve(device);

            // Initialize extra parameters.
            int numIterations = 0;
            int maxIterations = 20;
            double adjustmentFactor = 10;
            
            // Check if the initial solution produced a voltage at the interface that was too small.
            if (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage)
            {
                // The interace potential was too small. Iteratively increase the potential at the back of the device
                // until the interface voltage is too large (thus bracketing the correct value).
                while (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage)
                {
                    // Set the lower limit to the current potential at the back.
                    device.LowerLimit = device.Potential[0];

                    // Increase the potential at the back.
                    InitializeBoundary(device, adjustmentFactor * device.LowerLimit);

                    // Solve the device.
                    NumerovSolve(device);

                    // Increment the iteration counter.
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
                // We've bracketed the correct potential at the back. Store the upper limit. 
                device.UpperLimit = adjustmentFactor * device.Potential[0];
            }
            else
            {
                // The interace potential was too large. Iteratively decrease the potential at the back of the device
                // until the interface voltage is too small (thus bracketing the correct value).
                while (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge > device.Voltage)
                {
                    // Set the upper limit to the current potential at the back.
                    device.UpperLimit = device.Potential[0];

                    // Decrease the potential at the back.
                    InitializeBoundary(device, device.UpperLimit / adjustmentFactor);

                    // Solve the device.
                    NumerovSolve(device);

                    // Increment the iteration counter.
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
                // We've bracketed the correct potential at the back. Store the lower limit. 
                device.LowerLimit = device.Potential[0];
            }
        }

        /// <summary>
        /// Bracket the initial condition for the AC case.
        /// </summary>
        /// <param name="device">The device to be solved.</param>
        /// <param name="temperature">The temperature (in K).</param>
        /// <param name="frequency">The frequency (in Hz).</param>
        /// <param name="acVoltage">The AC voltage (in V).</param>
        internal static void Bracket(Device device, double temperature, double frequency, double acVoltage)
        {
            // Set the brackets to zero.
            device.UpperLimit = device.LowerLimit = 0;

            // Initialize the potenial and charge at the back of the device. Use the DC potential as a guess.
            InitializeBoundary(device, device.DCPotential[0]);

            // Solve for the AC potential and the AC charge density profiles.
            NoumerovSolve(device, temperature, frequency);

            // Initialize extra parameters.
            double adjustmentFactor = 10;

            // Check if the initial solution produced a voltage at the interface that was too small.
            if (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage - acVoltage)
            {
                // The interace potential was too small. Iteratively increase the potential at the back of the device
                // until the interface voltage is too large (thus bracketing the correct value).
                while (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge < device.Voltage - acVoltage)
                {
                    // Set the lower limit to the current potential at the back.
                    device.LowerLimit = device.Potential[0];

                    // Decrease the potential at the back.
                    InitializeBoundary(device, adjustmentFactor * device.LowerLimit);

                    // Solve the device.
                    NoumerovSolve(device, temperature, frequency);
                }
                // We've bracketed the correct potential at the back. Store the upper limit. 
                device.UpperLimit = adjustmentFactor * device.Potential[0];
            }
            else
            {
                // The interace potential was too large. Iteratively decrease the potential at the back of the device
                // until the interface voltage is too small (thus bracketing the correct value).
                while (device.Potential[device.Potential.Length - 1] / Constants.ElementaryCharge > device.Voltage - acVoltage)
                {
                    // Set the upper limit to the current potential at the back.
                    device.UpperLimit = device.Potential[0];

                    // Decrease the potential at the back.
                    InitializeBoundary(device, device.UpperLimit / adjustmentFactor);

                    // Solve the device.
                    NoumerovSolve(device, temperature, frequency);
                }
                // We've bracketed the correct potential at the back. Store the lower limit. 
                device.LowerLimit = device.Potential[0];
            }
        }

        /// <summary>
        /// Calculates the characteristic position x0 (in m) of the device.
        /// </summary>
        /// <param name="device">The device being solved.</param>
        /// <returns>Calculates the characteristic position x0 (in m) of the device.</returns>
        internal static double CalculateX0(Device device)
        {
            // Formula taken from Cohen and Lang's 1982 Phys. Rev. B paper.
            return Math.Sqrt(device.DielectricConstant / (Constants.ElementaryCharge * device.ShallowDopingDensity));
        }

        /// <summary>
        /// Finds DC solution to the Poisson's equation that matches the boundary conditions of the device.
        /// </summary>
        /// <param name="device">The device to be solved.</param>
        internal static void Solve(Device device)
        {
            // Initialize the interation counter and the tolerance.
            int numIterations = 0;
            double tolerance = 1e-6;

            device.Message += "Finding DC solution...\n";
            device.Message += "     Voltage    Target     Error      Tolerance\n";

            // Loop while the error in the interface voltage is greater than the tolerance.
            while (Math.Abs(device.Potential[device.NumberOfPositionPoints - 1] 
                / Constants.ElementaryCharge - device.Voltage) > tolerance)
            {
                // Set the potential at the back to half way between the upper and lower limit.
                InitializeBoundary(device, 0.5 * (device.UpperLimit + device.LowerLimit));

                // Solve Poisson's equation.
                NumerovSolve(device);

                device.Message += String.Format("{0,0:000}  {1,8:f7}  {2,8:f7}  {3,8:e2}  {4,8:e2}\n", 
                    numIterations,
                    device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge, 
                    device.Voltage,
                    Math.Abs(device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge - device.Voltage), 
                    tolerance);

                // Adjust the limits.
                if (device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge > device.Voltage)
                {
                    // The inteface voltage was too large. Reduce the upper limit to half way point.
                    device.UpperLimit = device.Potential[0];
                }
                else
                {
                    // The inteface voltage was too small. Increase the lower limit to half way point.
                    device.LowerLimit = device.Potential[0];
                }

                // Check if the solution is having trouble converging.
                if (numIterations >= 50 && device.NumRecursion < device.MaxRecursion)
                {
                    // If the thickness is too large the correct potential at the back cannot be represented by
                    // a double precision number (i.e., the double precision representation of the upper and lower 
                    // limit are identical). Alleviated by reducing the thickness.
                    device.NumRecursion++;
                    device.Message += "Maximum number of iterations reached!\nThickness is too large.\nReducing thickness by 20%.\n";
                    device.ChangeThickness(device.Thickness * 0.8);
                    Bracket(device);
                    Solve(device); // Recursive call.
                    break;
                }

                // Increment the interation counter.
                numIterations++;
            }

            // Check if the device is too thin.
            if (device.ChargeDensity[0] / device.ChargeDensity.Max() > 1e-4)
            {
                // If the device is too thin, then there will be no neutral region at the back of the
                // device. This violates the assumptions of the algorithm, thus the thickness must be increased.
                // A charge density 10,000x smaller than the maximum value is considered sufficiently neutral.
                device.NumRecursion++;
                device.Message += "Thickness is too small.\nIncreasing thickness by 10%.\n";
                device.ChangeThickness(device.Thickness * 1.1);
                Bracket(device);
                Solve(device); // Recursive call.
            }
        }

        /// <summary>
        /// Finds AC solution to the Poisson's equation that matches the boundary conditions of the device.
        /// </summary>
        /// <param name="device">The device to be solved.</param>
        /// <param name="temperature">The temperature (in K) of the device.</param>
        /// <param name="frequency">The frequency (in Hz) of the AC perturbation.</param>
        /// <param name="acVoltage">The AC perturbation (in V) applied to the device. Forward bias is positive.</param>
        internal static void Solve(Device device, double temperature, double frequency, double acVoltage)
        {
            // Initialize the interation counter and the tolerance.
            int numIterations = 0;
            double tolerance = 1e-6;

            // Loop while the error in the interface voltage is greater than the tolerance.
            while (Math.Abs(device.Potential[device.NumberOfPositionPoints - 1] 
                / Constants.ElementaryCharge - (device.Voltage - acVoltage)) > tolerance)
            {
                // Set the potential at the back to half way between the upper and lower limit.
                InitializeBoundary(device, 0.5 * (device.UpperLimit + device.LowerLimit));

                // Solve Poisson's equation.
                NoumerovSolve(device, temperature, frequency);

                // Adjust the limits.
                if (device.Potential[device.NumberOfPositionPoints - 1] / Constants.ElementaryCharge > device.Voltage - acVoltage)
                {
                    // The inteface voltage was too large. Reduce the upper limit to half way point.
                    device.UpperLimit = device.Potential[0];
                }
                else
                {
                    // The inteface voltage was too small. Increase the lower limit to half way point.
                    device.LowerLimit = device.Potential[0];
                }

                // Check if the solution is having trouble converging.
                if (numIterations >= 50 && device.NumRecursion < device.MaxRecursion)
                {
                    // If the thickness is too large the correct potential at the back cannot be represented by
                    // a double precision number (i.e., the double precision representation of the upper and lower 
                    // limit are identical). Alleviated by reducing the thickness.
                    device.NumRecursion++;
                    device.ChangeThickness(device.Thickness * 0.8);
                    Bracket(device);
                    Solve(device);
                    device.CopyResultsToDCArrays();
                    Bracket(device, temperature, frequency, acVoltage);
                    Solve(device, temperature, frequency, acVoltage); // Recursive call.
                    break;
                }

                // Increment the interation counter.
                numIterations++;
            }

            // Check if the device is too thin.
            if (device.ChargeDensity[0] / device.ChargeDensity.Max() > 1e-4)
            {
                // If the device is too thin, then there will be no neutral region at the back of the
                // device. This violates the assumptions of the algorithm, thus the thickness must be increased.
                // A charge density 10,000x smaller than the maximum value is considered sufficiently neutral.
                device.NumRecursion++;
                device.ChangeThickness(device.Thickness * 1.1);
                Bracket(device);
                Solve(device);
                device.CopyResultsToDCArrays();
                Bracket(device, temperature, frequency, acVoltage);
                Solve(device, temperature, frequency, acVoltage);
            }
        }

        /// <summary>
        /// Calculates the demarcation energy of the device.
        /// </summary>
        /// <param name="device">The device in question.</param>
        /// <param name="temperature">The temperature (in K) of the device.</param>
        /// <param name="frequency">The frequency (in Hz) of the AC perturbation.</param>
        /// <returns>The demarcation energy (in J).</returns>
        internal static double CalculateDemarcationEnergy(Device device, double temperature, double frequency)
        {
            return Constants.BoltzmannConstant * temperature * Math.Log(device.ThermalEmissionPrefactor * temperature * temperature / frequency);
        }

        #endregion [Internal Methods]
    }
}
