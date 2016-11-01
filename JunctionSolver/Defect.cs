using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    public class Defect
    {
        public string Label { get; set; }
        public string PositionDependence { get; set; }

        public double Energy { get; set; }
        public double FWHM { get; set; }
        public double Magnitude { get; set; }

        public double[][] DensityOfStates { get; set; }

        public Defect(Device device, double energy, double fwhm, double magnitude, string label, string positionDependence)
        {
            Energy = energy;
            FWHM = fwhm;
            Magnitude = magnitude;
            Label = label;
            PositionDependence = positionDependence;

            DensityOfStates = new double[device.NumberOfPositionPoints][];
            
            if (positionDependence == "Constant" || positionDependence == "constant" || positionDependence == "1") 
            {
                for (int i = 0; i < device.NumberOfPositionPoints; i++)
                {
                    DensityOfStates[i] = Utils.Gaussian(device.Energy, Energy, FWHM / 2.3548);

                    for (int j = 0; j < device.NumberOfEnergyPoints; j++)
                    {
                        DensityOfStates[i][j] *= Magnitude;
                    }
                }
            }
            else
            {
                NCalc.Expression posDep = new NCalc.Expression(PositionDependence);

                for (int i = 0; i < device.NumberOfPositionPoints; i++)
                {
                    DensityOfStates[i] = Utils.Gaussian(device.Energy, Energy, FWHM / 2.3548);

                    posDep.Parameters["x"] = device.FlippedPosition[i];

                    double posDepFactor = Math.Max(0.0, Convert.ToDouble(posDep.Evaluate()));

                    for (int j = 0; j < device.NumberOfEnergyPoints; j++)
                    {
                        DensityOfStates[i][j] *= Magnitude * posDepFactor;
                    }
                }
            }
        }

        /*
        public Defect(Device device, double energy, double fwhm, double magnitude, string label)
        {
            Energy = energy;
            FWHM = fwhm;
            Magnitude = magnitude;
            Label = label;
            PositionDependence = "Constant";

            DensityOfStates = new double[device.NumberOfPositionPoints][];

            for (int i = 0; i < device.NumberOfPositionPoints; i++)
            {
                DensityOfStates[i] = Utils.Gaussian(device.Energy, Energy, FWHM / 2.3548);

                for (int j = 0; j < device.NumberOfEnergyPoints; j++)
                {
                    DensityOfStates[i][j] *= Magnitude;
                }
            }
        }
        */
    }
}
