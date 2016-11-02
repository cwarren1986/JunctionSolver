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
        public double[][] ChargeDensityTable { get; set; }

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
        /// <param name="numberOfPositionPoints">The number of points in the position grid.</param>
        /// <param name="numberOfEnergyPoints">The number of points in the energy grid.</param>
        /// <param name="dielectricConstant">The dielectric constant of the device in F/m.</param>
        /// <param name="thickness">The thickness of the model in m.</param>
        /// <param name="bandGap">The bandgap of the device in J.</param>
        /// <param name="neutralRegionFermiLevel">The neutral region Fermi level in J. Should be negative.</param>
        /// <param name="shallowDopingDensity">The shallow doping density of the device in /m^3.</param>
        /// <param name="thermalEmissionPrefactor">The thermal emission prefactor of the device in Hz/K^2.</param>
        /// <param name="builtInVoltage">The built-in voltage of the device in V.</param>
        /// <param name="appliedVoltage">The DC voltage applied to the device in V.</param>
        /// <param name="defectList">A list containing Defect objects to be included in the model.</param>
        public Device(
            int numberOfPositionPoints,
            int numberOfEnergyPoints,
            double dielectricConstant,
            double thickness,
            double bandGap,
            double neutralRegionFermiLevel,
            double shallowDopingDensity,
            double thermalEmissionPrefactor,
            double builtInVoltage,
            double appliedVoltage,
            List<Defect> defectList
            )
        {
            // Copy input parameters to the appropriate properties and fields.
            NumberOfPositionPoints = numberOfPositionPoints;
            NumberOfEnergyPoints = numberOfEnergyPoints;
            DielectricConstant = dielectricConstant;
            Thickness = thickness;
            BandGap = bandGap;
            NeutralRegionFermiLevel = neutralRegionFermiLevel;
            ShallowDopingDensity = shallowDopingDensity;
            ThermalEmissionPrefactor = thermalEmissionPrefactor;
            BuiltInVoltage = builtInVoltage;
            AppliedVoltage = appliedVoltage;

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
            Potential = new double[NumberOfPositionPoints];
            ChargeDensity = new double[NumberOfPositionPoints];
            FermiLevel = new double[NumberOfPositionPoints];

            DCPotential = new double[NumberOfPositionPoints];
            DCChargeDensity = new double[NumberOfPositionPoints];
            DCFermiLevel = new double[NumberOfPositionPoints];

            Energy = new double[NumberOfEnergyPoints];
            DensityOfStates = new double[NumberOfPositionPoints][];
            ChargeDensityTable = new double[NumberOfPositionPoints][];

            // Fill position space arrays, define second dimension of 2D arrays.
            for (int i = 0; i < NumberOfPositionPoints; i++)
            {
                Position[i] = PositionSpacing * i;
                FlippedPosition[NumberOfPositionPoints - 1 - i] = Position[i];
                FermiLevel[i] = NeutralRegionFermiLevel;
                DensityOfStates[i] = new double[NumberOfEnergyPoints];
                ChargeDensityTable[i] = new double[NumberOfEnergyPoints];
            }
            
            // Fill energy space arrays.
            for (int i = 0; i < NumberOfEnergyPoints; i++)
            {
                Energy[i] = NeutralRegionFermiLevel + i * EnergySpacing;
            }

            // Initialize the DefectList.
            DefectList = new List<Defect>();

            // Fill DefectList.
            for (int i = 0; i < defectList.Count; i++)
            {
                // If the defect is the shallow dopant, explicitly force it to be
                // centered on the Fermi level, and to have the correct density.
                if (defectList[i].Label == "Shallow Dopant")
                {
                    DefectList.Add(new Defect(
                        this,
                        NeutralRegionFermiLevel,
                        defectList[i].FWHM,
                        2 * ShallowDopingDensity * Constants.ElementaryCharge, // Factor of 2 because only half the defect is integrated over.
                        defectList[i].Label,
                        defectList[i].PositionDependence
                    ));
                }
                else
                {
                    DefectList.Add(new Defect(
                        this,
                        defectList[i].Energy,
                        defectList[i].FWHM,
                        defectList[i].Magnitude,
                        defectList[i].Label,
                        defectList[i].PositionDependence
                    ));
                }
            }

            // If the defect list is empty add the shallow dopant.
            if (DefectList == null || DefectList.Count == 0)
            {
                DefectList.Add(new Defect(
                this,
                NeutralRegionFermiLevel,
                30e-3 * Constants.ElementaryCharge,
                2 * ShallowDopingDensity * Constants.ElementaryCharge,
                "Shallow Dopant",
                "Constant"
                ));
            }

            // Calculate the density of states.
            Utils.CalculateDensityOfStates(this);

            // Precalculate the charge density integrals.
            Utils.RhoTable(this);

            // Calculate the characteristic length.
            X0 = Utils.Calcx0(this);
        }

        #endregion [Constructor]

        #region [Public Methods]

        /// <summary>
        /// Copies model results to the DC arrays.
        /// </summary>
        public void CopyResultsToDCArrays()
        {
            Potential.CopyTo(DCPotential, 0);
            ChargeDensity.CopyTo(DCChargeDensity, 0);
            FermiLevel.CopyTo(DCFermiLevel, 0);
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

            // Recalculate the density of states and the charge density integral table.
            Utils.CalculateDensityOfStates(this);

            Utils.RhoTable(this);

            // Recalculate the characteristic length.
            X0 = Utils.Calcx0(this);
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

        #endregion [Public Methods]
    }
}
