using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    /// <summary>
    /// A class for holding the properties of a device.
    /// </summary>
    public class DeviceParameters
    {

        #region [Public Properties]

        /// <summary>
        /// The number of points in the position grid.
        /// </summary>
        public int NumberOfPositionPoints { get; set; }

        /// <summary>
        /// The number of points in the energy grid.
        /// </summary>
        public int NumberOfEnergyPoints { get; set; }

        /// <summary>
        /// The dielectric constant (in F/m) of the device.
        /// </summary>
        public double DielectricConstant { get; set; }

        /// <summary>
        /// The thickness (in m) of the device model.
        /// </summary>
        public double Thickness { get; set; }

        /// <summary>
        /// The bandgap (in J) of the device.
        /// </summary>
        public double BandGap { get; set; }

        /// <summary>
        /// The Fermi level (in J) in the neutral bulk.
        /// </summary>
        /// <remarks>Should be negative since the potential is zero in the neutral region.</remarks>
        public double NeutralRegionFermiLevel { get; set; }

        /// <summary>
        /// The shallow doping density (in /m^3) of the device.
        /// </summary>
        public double ShallowDopingDensity { get; set; }

        /// <summary>
        /// The thermal emmission prefactor (in Hz/K^2) of the device.
        /// </summary>
        public double ThermalEmissionPrefactor { get; set; }

        /// <summary>
        /// The built-in voltage (in V) of the device.
        /// </summary>
        /// <remarks>Forward bias is positive.</remarks>
        public double BuiltInVoltage { get; set; }

        /// <summary>
        /// The DC voltage (in V) applied to the device.
        /// </summary>
        /// <remarks>Forward bias is positive.</remarks>
        public double AppliedVoltage { get; set; }

        /// <summary>
        /// A list of DefectParameters objects that define the density of states of the device.
        /// </summary>
        public List<DefectParameters> DefectParameterList { get; set; }

        #endregion [Public Properties]

        #region [Constructor]

        /// <summary>
        /// The constructor.
        /// </summary>
        /// <param name="numberOfPositionPoints">The number of points in the position grid.</param>
        /// <param name="numberOfEnergyPoints">The number of points in the energy grid.</param>
        /// <param name="dielectricConstant">The dielectric constant (in F/m) of the device.</param>
        /// <param name="thickness">The thickness (in m) of the device model.</param>
        /// <param name="bandGap">The bandgap (in J) of the device.</param>
        /// <param name="neutralRegionFermiLevel">The Fermi level (in J) in the neutral bulk.</param>
        /// <param name="shallowDopingDensity">The shallow doping density (in /m^3) of the device.</param>
        /// <param name="thermalEmissionPrefactor">The thermal emmission prefactor (in Hz/K^2) of the device.</param>
        /// <param name="builtInVoltage">The built-in voltage (in V) of the device.</param>
        /// <param name="appliedVoltage">The DC voltage (in V) applied to the device.</param>
        /// <param name="defectParameterList">A list of DefectParameters objects that define the density of states of the device.</param>
        public DeviceParameters(
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
            List<DefectParameters> defectParameterList
            )
        {
            // Copy the parameters to their corresponding properties.
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
            DefectParameterList = defectParameterList;
        }

        #endregion [Constructor]
    }
}
