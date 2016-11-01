using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    public class Device
    {
        public string Message { get; set; }

        public int NumberOfPositionPoints { get; set; }
        public int NumberOfEnergyPoints { get; set; }

        public double DielectricConstant { get; set; }
        public double Thickness { get; set; }
        public double BandGap { get; set; }
        public double NeutralRegionFermiLevel { get; set; }
        public double BuiltInVoltage { get; set; }
        public double AppliedVoltage { get; set; }
        public double ShallowDopingDensity { get; set; }
        public double ThermalEmissionPrefactor { get; set; }

        public double[] Position { get; set; }
        public double[] FlippedPosition { get; set; }
        public double[] Potential { get; set; }
        public double[] ChargeDensity { get; set; }
        public double[] FermiLevel { get; set; }
        public double[] Energy { get; set; }
        public double[][] DensityOfStates { get; set; }
        public double[][] ChargeDensityTable { get; set; }

        public double[] DCPotential { get; set; }
        public double[] DCChargeDensity { get; set; }
        public double[] DCFermiLevel { get; set; }

        public List<Defect> DefectList { get; set; }

        internal double UpperLimit;
        internal double LowerLimit;
        internal double Voltage;
        internal double PositionSpacing;
        internal double PositionSpacingSquared;
        internal double HalfBandGap;
        internal double EnergySpacing;
        internal double X0;
        internal double ElectricFluxConstant;
        internal int NumRecursion = 0;
        internal int MaxRecursion = 10;

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
            NumberOfPositionPoints = numberOfPositionPoints;
            NumberOfEnergyPoints = numberOfEnergyPoints;
            DielectricConstant = dielectricConstant;
            Thickness = thickness;
            BandGap = bandGap;
            HalfBandGap = 0.5 * BandGap;
            NeutralRegionFermiLevel = neutralRegionFermiLevel;
            ShallowDopingDensity = shallowDopingDensity;
            ThermalEmissionPrefactor = thermalEmissionPrefactor;
            BuiltInVoltage = builtInVoltage;
            AppliedVoltage = appliedVoltage;
            Voltage = BuiltInVoltage - AppliedVoltage;
            ElectricFluxConstant = Constants.ElementaryCharge / DielectricConstant;

            EnergySpacing = -HalfBandGap / (NumberOfEnergyPoints - 1);
            PositionSpacing = Thickness / (NumberOfPositionPoints - 1);
            PositionSpacingSquared = PositionSpacing * PositionSpacing;

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

            for (int i = 0; i < NumberOfPositionPoints; i++)
            {
                Position[i] = PositionSpacing * i;
                FlippedPosition[NumberOfPositionPoints - 1 - i] = Position[i];
                FermiLevel[i] = NeutralRegionFermiLevel;
                DensityOfStates[i] = new double[NumberOfEnergyPoints];
                ChargeDensityTable[i] = new double[NumberOfEnergyPoints];
            }
            
            for (int i = 0; i < NumberOfEnergyPoints; i++)
            {
                Energy[i] = NeutralRegionFermiLevel + i * EnergySpacing;
            }

            DefectList = new List<Defect>();//defectList;
            for (int i = 0; i < defectList.Count; i++)
            {
                if (defectList[i].Label == "Shallow Dopant")
                {
                    DefectList.Add(new Defect(
                        this,
                        NeutralRegionFermiLevel,
                        defectList[i].FWHM,
                        2 * ShallowDopingDensity * Constants.ElementaryCharge,
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

            Utils.Calcg(this);

            Utils.RhoTable(this);

            X0 = Utils.Calcx0(this);
        }

        public void CopyResultsToDCArrays()
        {
            Potential.CopyTo(DCPotential, 0);
            ChargeDensity.CopyTo(DCChargeDensity, 0);
            FermiLevel.CopyTo(DCFermiLevel, 0);
        }

        public void ChangeAppliedVoltage(double appliedVoltage)
        {
            AppliedVoltage = appliedVoltage;
            Voltage = BuiltInVoltage - AppliedVoltage;
        }

        public void ChangeThickness(double thickness)
        {
            Thickness = thickness;

            PositionSpacing = Thickness / (NumberOfPositionPoints - 1);

            PositionSpacingSquared = PositionSpacing * PositionSpacing;

            for (int i = 0; i < NumberOfPositionPoints; i++)
            {
                Position[i] = PositionSpacing * i;
                FlippedPosition[NumberOfPositionPoints - 1 - i] = Position[i];
            }

            Utils.Calcg(this);

            Utils.RhoTable(this);

            X0 = Utils.Calcx0(this);
        }

        public void Solve()
        {
            NumRecursion = 0;

            Utils.Bracket(this);

            Message += "Solution bracketed.\n";

            NumRecursion = 0;

            Utils.Solve(this);

            Message += "Acceptable solution found.\n";

            CopyResultsToDCArrays();
        }

        public void Solve(double temperature, double frequency, double acVoltage)
        {
            NumRecursion = 0;

            Utils.Bracket(this, temperature, frequency, acVoltage);

            NumRecursion = 0;

            Utils.Solve(this, temperature, frequency, acVoltage);
        }
    }
}
