using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    // A class that holds and manages the properties of a device (i.e. a one-sided p-n junction).  
    public class Device
    {
        #region [Public Propeties]

        /// <summary>
        /// Holds messages, warnings, etc. so they can be printed if needed. 
        /// </summary>
        public string Message { get; set; }

        /// <summary>
        /// The number of points in position space used by the model.
        /// </summary>
        public int NumberOfPositionPoints { get; set; }

        /// <summary>
        /// The number of points in energy space used by the model.
        /// </summary>
        public int NumberOfEnergyPoints { get; set; }

        /// <summary>
        /// The dielectric constant (in F/m) of the lightly doped side of the junction.
        /// </summary>
        public double DielectricConstant { get; set; }

        /// <summary>
        /// The total thickness (in m) of the model.
        /// </summary>
        public double Thickness { get; set; }

        /// <summary>
        /// The bandgap (in J) of the device.
        /// </summary>
        public double BandGap { get; set; }

        /// <summary>
        /// The fermi level (in J) in the neutral region of the device.
        /// Since the pontential is defined to be zero in the neutral bulk,
        /// and the built-in voltage is defined to be positve, this should 
        /// be a negative number.
        /// </summary>
        public double NeutralRegionFermiLevel { get; set; }

        /// <summary>
        /// The built in voltage (in V) of the device. Defined to be positive.
        /// </summary>
        public double BuiltInVoltage { get; set; }

        /// <summary>
        /// The DC voltage (in V) applied to the device. Forward bias is positive.
        /// </summary>
        public double AppliedVoltage { get; set; }

        /// <summary>
        /// The shallow doping density (in /m^3) of the device.
        /// </summary>
        public double ShallowDopingDensity { get; set; }

        /// <summary>
        /// The thermal emission prefactor (in Hz/K^2) of the device.
        /// </summary>
        public double ThermalEmissionPrefactor { get; set; }

        /// <summary>
        /// The position grid (in m) of the model.
        /// The position is zero at the back of the device.
        /// </summary>
        public double[] Position { get; set; }

        /// <summary>
        /// The position grid (in m) of the model subtracted from the thickness
        /// such that the position is zero at the interface instead of the back.
        /// </summary>
        public double[] FlippedPosition { get; set; }

        /// <summary>
        /// Array for storing the potential profile (in J). Defined to be zero in the neutral bulk.
        /// Should be strictly positive.
        /// </summary>
        public double[] Potential { get; set; }

        /// <summary>
        /// Array for storing the charge density profile (in C/m^3). Defined to be zero in the bulk.
        /// Should be strictly positive.
        /// </summary>
        public double[] ChargeDensity { get; set; }

        /// <summary>
        /// Array for storing the fermi level profile (in J).
        /// </summary>
        public double[] FermiLevel { get; set; }

        /// <summary>
        /// The energy grid (in J) of the model.
        /// </summary>
        public double[] Energy { get; set; }

        /// <summary>
        /// Array for storing the density of states (in /m^3eV) of the model as a function of
        /// position (first index) and energy (second index).
        /// </summary>
        public double[][] DensityOfStates { get; set; }

        /// <summary>
        /// Array for storing precalculated charge density integrals (in C/m^3) as a function of
        /// position (first index) and energy (second index).
        /// </summary>
        public double[][] ChargeDensityIntegralTable { get; set; }

        /// <summary>
        /// Array for storing the DC potential profile (in J).
        /// </summary>
        public double[] DCPotential { get; set; }

        /// <summary>
        /// Array for storing the DC charge density profile (in C/m^3).
        /// </summary>
        public double[] DCChargeDensity { get; set; }

        /// <summary>
        /// Array for storing the DC fermi level profile (in J).
        /// </summary>
        public double[] DCFermiLevel { get; set; }

        public List<DefectParameters> DefectParameterList { get; set; }

        /// <summary>
        /// List for holding Defect objects that are used to construct the
        /// density of states.
        /// </summary>
        public List<Defect> DefectList { get; set; }

        #endregion [Public Properties]

        #region [Internal Fields]

        /// <summary>
        /// A value of the pontential (in J) at the back 
        /// of the device which is larger than the correct value.
        /// Used to bracket the correct solution.
        /// </summary>
        internal double UpperLimit;

        /// <summary>
        /// A value of the pontential (in J) at the back 
        /// of the device which is smaller than the correct value.
        /// Used to bracket the correct solution.
        /// </summary>
        internal double LowerLimit;

        /// <summary>
        /// The target voltage (in V) at the interface including both the
        /// built-in voltage and the DC applied voltage. Does not include the
        /// applied AC voltage.
        /// </summary>
        internal double Voltage;

        /// <summary>
        /// The spacing (in m) of the position grid.
        /// </summary>
        internal double PositionSpacing;

        /// <summary>
        /// The squared spacing (in m^2) of the position grid.
        /// </summary>
        internal double PositionSpacingSquared;

        /// <summary>
        /// The position spacing (in m) of the position grid for the DC solution.
        /// </summary>
        internal double DCPositionSpacing;

        /// <summary>
        /// One half of the bandgap (in J).
        /// </summary>
        internal double HalfBandGap;

        /// <summary>
        /// The spacing of the energy grid (in J).
        /// </summary>
        internal double EnergySpacing;

        /// <summary>
        /// The characteristic length of the potential profile for the
        /// case when the density of states in constant in both position
        /// and energy space. Can be determined analytically. Useful for
        /// guessing the potential at the back of the device.
        /// </summary>
        internal double X0;

        /// <summary>
        /// The electric flux constant (in Vm) of the deivce (saves repeatedly calculating q/ε).
        /// </summary>
        internal double ElectricFluxConstant;

        /// <summary>
        /// The current recursion depth of the solver.
        /// </summary>
        internal int NumRecursion = 0;

        /// <summary>
        /// The maximum recursion depth of the solver.
        /// </summary>
        internal int MaxRecursion = 10;

        #endregion [Internal Fields]

        #region [Constructor]

        /// <summary>
        /// The constructor.
        /// </summary>
        /// <param name="deviceParameters">Object that stores the device parameters</param>
        public Device(DeviceParameters deviceParameters)
        {
            // Copy input parameters to the appropriate properties and fields.
            NumberOfPositionPoints = deviceParameters.NumberOfPositionPoints;
            NumberOfEnergyPoints = deviceParameters.NumberOfEnergyPoints;
            DielectricConstant = deviceParameters.DielectricConstant;
            Thickness = deviceParameters.Thickness;
            BandGap = deviceParameters.BandGap;
            NeutralRegionFermiLevel = deviceParameters.NeutralRegionFermiLevel;
            ShallowDopingDensity = deviceParameters.ShallowDopingDensity;
            ThermalEmissionPrefactor = deviceParameters.ThermalEmissionPrefactor;
            BuiltInVoltage = deviceParameters.BuiltInVoltage;
            AppliedVoltage = deviceParameters.AppliedVoltage;
            DefectParameterList = deviceParameters.DefectParameterList;

            // Calculate other properties and fields.
            HalfBandGap = 0.5 * BandGap;
            Voltage = BuiltInVoltage - AppliedVoltage;
            ElectricFluxConstant = Constants.ElementaryCharge / DielectricConstant;
            EnergySpacing = -HalfBandGap / (NumberOfEnergyPoints - 1);
            PositionSpacing = Thickness / (NumberOfPositionPoints - 1);
            PositionSpacingSquared = PositionSpacing * PositionSpacing;

            // Initialize arrays.
            Position = new double[NumberOfPositionPoints];
            FlippedPosition = new double[NumberOfPositionPoints];
            //DCFlippedPosition = new double[NumberOfPositionPoints];
            Potential = new double[NumberOfPositionPoints];
            ChargeDensity = new double[NumberOfPositionPoints];
            FermiLevel = new double[NumberOfPositionPoints];

            //DCPotential = new double[NumberOfPositionPoints];
            //DCChargeDensity = new double[NumberOfPositionPoints];
            //DCFermiLevel = new double[NumberOfPositionPoints];

            Energy = new double[NumberOfEnergyPoints];
            DensityOfStates = new double[NumberOfPositionPoints][];
            ChargeDensityIntegralTable = new double[NumberOfPositionPoints][];

            // Fill position space arrays, define second dimension of 2D arrays.
            for (int i = 0; i < NumberOfPositionPoints; i++)
            {
                Position[i] = PositionSpacing * i;
                FlippedPosition[NumberOfPositionPoints - 1 - i] = Position[i];
                FermiLevel[i] = NeutralRegionFermiLevel;
                DensityOfStates[i] = new double[NumberOfEnergyPoints];
                ChargeDensityIntegralTable[i] = new double[NumberOfEnergyPoints];
            }
            
            // Fill energy space arrays.
            for (int i = 0; i < NumberOfEnergyPoints; i++)
            {
                Energy[i] = NeutralRegionFermiLevel + i * EnergySpacing;
            }

            // Initialize the DefectList.
            DefectList = new List<Defect>();
            
            // Fill DefectList.
            for (int i = 0; i < DefectParameterList.Count; i++)
            {
                DefectList.Add(new Defect(this, DefectParameterList[i]));
            }

            // Calculate the density of states.
            Utils.CalculateDensityOfStates(this);

            // Precalculate the charge density integrals.
            Utils.CalculateChargeDensityIntegralTable(this);

            // Calculate the characteristic length.
            X0 = Utils.CalculateX0(this);
        }

        #endregion [Constructor]

        #region [Public Methods]

        /// <summary>
        /// Copies model results to the DC arrays.
        /// </summary>
        public void CopyResultsToDCArrays()
        {
            DCPotential = new double[NumberOfPositionPoints];
            DCChargeDensity = new double[NumberOfPositionPoints];
            DCFermiLevel = new double[NumberOfPositionPoints];

            Potential.CopyTo(DCPotential, 0);
            ChargeDensity.CopyTo(DCChargeDensity, 0);
            FermiLevel.CopyTo(DCFermiLevel, 0);
            DCPositionSpacing = PositionSpacing;
        }

        /// <summary>
        /// Changes the voltage applied to the device.
        /// Recalculates other parameters as required.
        /// </summary>
        /// <param name="appliedVoltage">The voltage (in V) applied to the device.</param>
        public void ChangeAppliedVoltage(double appliedVoltage)
        {
            // Update applied voltage.
            AppliedVoltage = appliedVoltage;

            // Recalculate Voltage property.
            Voltage = BuiltInVoltage - AppliedVoltage;
        }

        /// <summary>
        /// Changes the thickness (in m) of the model.
        /// Recalculates other parameters as requried.
        /// </summary>
        /// <param name="thickness">The thickness (in m) of the model.</param>
        public void ChangeThickness(double thickness)
        {
            // Update thickness.
            Thickness = thickness;

            // Recalculate position spacing and position spacing squared.
            PositionSpacing = Thickness / (NumberOfPositionPoints - 1);

            PositionSpacingSquared = PositionSpacing * PositionSpacing;

            // Refill the position and flipped position grids.
            for (int i = 0; i < NumberOfPositionPoints; i++)
            {
                Position[i] = PositionSpacing * i;
                FlippedPosition[NumberOfPositionPoints - 1 - i] = Position[i];
            }

            // Initialize the DefectList.
            DefectList = new List<Defect>();

            // Fill DefectList.
            for (int i = 0; i < DefectParameterList.Count; i++)
            {
                DefectList.Add(new Defect(this, DefectParameterList[i]));
            }

            // Recalculate the density of states and the charge density integral table.
            Utils.CalculateDensityOfStates(this);

            Utils.CalculateChargeDensityIntegralTable(this);

            // Recalculate the characteristic length.
            X0 = Utils.CalculateX0(this);
        }

        /// <summary>
        /// Calculates the DC potential that solves Poisson's equation in the device.
        /// </summary>
        public void Solve()
        {
            // Reset the recursion counter.
            NumRecursion = 0;

            // Bracket the solution.
            Utils.Bracket(this);

            Message += "Solution bracketed.\n";

            // Reset the recursion counter.
            NumRecursion = 0;

            // Find the solution.
            Utils.Solve(this);

            Message += "Acceptable solution found.\n";

            // Store results in DC arrays.
            CopyResultsToDCArrays();
        }

        /// <summary>
        /// Calculates the AC potential that solves Poisson's equation in the device. 
        /// </summary>
        /// <param name="temperature">The temperature (in K) of the deivce.</param>
        /// <param name="frequency">The frequency (in Hz) of the AC voltage applied to the device.</param>
        /// <param name="acVoltage">The AC voltage (in V) applied to the device. Forward bias is positive.</param>
        public void Solve(double temperature, double frequency, double acVoltage)
        {
            // Reset the recursion counter.
            NumRecursion = 0;

            // Bracket the solution.
            Utils.Bracket(this, temperature, frequency, acVoltage);

            // Reset the recursion counter.
            NumRecursion = 0;

            // Find the solution.
            Utils.Solve(this, temperature, frequency, acVoltage);
        }

        /// <summary>
        /// Calculates the capacitance (in F/m^2) of the device.
        /// </summary>
        /// <returns></returns>
        public double CalculateCapacitance()
        {
            // Calculate the change in voltage.
            double dV = (Potential[NumberOfPositionPoints - 1]
                - DCPotential[NumberOfPositionPoints - 1])
                / Constants.ElementaryCharge;

            // Calculate the charge in the AC and DC cases.
            double totalChargeAC = 0;
            double totalChargeDC = 0;
            for (int i = 0; i < NumberOfPositionPoints; i++)
            {
                totalChargeAC += ChargeDensity[i];
                totalChargeDC += DCChargeDensity[i];
            }

            // Calcualate the change in charge.
            double dQ = totalChargeAC * PositionSpacing - totalChargeDC * DCPositionSpacing;

            // Return the capacitance.
            return dQ / dV;
        }

        /// <summary>
        /// Calcualtes the demarcation energy (in J) of the device.
        /// </summary>
        /// <param name="temperature">The temperature (in K).</param>
        /// <param name="frequency">The frequency (in Hz).</param>
        /// <returns>The demarcation energy (in J).</returns>
        public double CalculateDemarcationEnergy(double temperature, double frequency)
        {
            return Utils.CalculateDemarcationEnergy(this, temperature, frequency);
        }

        #endregion [Public Methods]
    }
}
