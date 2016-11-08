using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    /// <summary>
    /// A class for holding defect properties.
    /// </summary>
    public class Defect
    {
        #region [Private Fields]

        private double fwhmPerStandardDeviation;

        #endregion [Private Fields]

        #region [Public Properties]

        /// <summary>
        /// The parameters of the defect.
        /// </summary>
        public DefectParameters DefectParameters { get; set; }

        /// <summary>
        /// The density of states of the defect as a function of position (1st index) and energy (2nd index).
        /// </summary>
        /// <remarks>The position and energy grid is defined by the device.</remarks>
        public double[][] DensityOfStates { get; set; }

        #endregion [Public Properties]

        #region [Constructor]

        /// <summary>
        /// The constructor.
        /// </summary>
        /// <param name="device">The device being solved.</param>
        /// <param name="defectParameters">The parameters of the defect.</param>
        public Defect(Device device, DefectParameters defectParameters)
        {
            // Define the DefectParameters property.
            DefectParameters = defectParameters;

            // Define factor that converts FWHM to standard deviation.
            fwhmPerStandardDeviation = 2 * Math.Sqrt(2 * Math.Log(2));

            // Initialize the 1st index of the density of states.
            DensityOfStates = new double[device.NumberOfPositionPoints][];
            
            // Fill the density of states array.
            if (DefectParameters.PositionDependence == "Constant" 
                || DefectParameters.PositionDependence == "constant" 
                || DefectParameters.PositionDependence == "1") 
            {
                // The position dependence is constant, allowing us to simplify filling the array.
                // Namely, we can precalculate the energy dependence of the DOS since it'll be the same
                // at every position.
                double[] dos = Utils.Gaussian(device.Energy, DefectParameters.Energy, 
                    DefectParameters.FWHM / fwhmPerStandardDeviation, DefectParameters.Magnitude);

                // Loop over position grid.
                for (int i = 0; i < device.NumberOfPositionPoints; i++)
                {
                    // Fill a row of DensityOfStates.
                    DensityOfStates[i] = dos;
                }
            }
            else
            {
                // The position dependence is not constant. Evaluate the position dependence expression.
                NCalc.Expression posDep = new NCalc.Expression(DefectParameters.PositionDependence);

                // Loop over the position grid.
                for (int i = 0; i < device.NumberOfPositionPoints; i++)
                {
                    // Set the position variable in the position dependence to the flipped position value.
                    // The user only ever sees the flipped position, so they define things in terms of it.
                    posDep.Parameters["x"] = device.FlippedPosition[i];

                    // Calculate the position dependence factor, while not allowing it to be negative.
                    double posDepFactor = Math.Max(0.0, Convert.ToDouble(posDep.Evaluate()));

                    // Fill a row of the density of states.
                    DensityOfStates[i] = Utils.Gaussian(device.Energy, DefectParameters.Energy, 
                        DefectParameters.FWHM / fwhmPerStandardDeviation, DefectParameters.Magnitude * posDepFactor);
                }
            }
        }

        #endregion [Constructor]
    }
}
