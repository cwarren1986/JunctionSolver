using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    /// <summary>
    /// A class for holding the parameters of a Gaussian defect.
    /// </summary>
    public class DefectParameters
    {

        #region [Public Properties]

        /// <summary>
        /// A label for the defect.
        /// </summary>
        public string Label { get; set; }

        /// <summary>
        /// The position dependence of the defect.
        /// </summary>
        public string PositionDependence { get; set; }

        /// <summary>
        /// The energy (in J) at which the defect is centered.
        /// </summary>
        public double Energy { get; set; }

        /// <summary>
        /// The full width a half maximum (in J) of the defect.
        /// </summary>
        public double FWHM { get; set; }

        /// <summary>
        /// The magnitude (in /m^3) of the defect, i.e., the density you get by integrating over energy.
        /// </summary>
        public double Magnitude { get; set; }

        #endregion [Public Properties]

        #region [Constructor]

        /// <summary>
        /// The constructor.
        /// </summary>
        /// <param name="energy">The energy of the defect (in J).</param>
        /// <param name="fwhm">The full width at half maximum (in J) of the defect.</param>
        /// <param name="magnitude">The magnitude of the defect (in /m^3) of the defect.</param>
        /// <param name="label">The label of the defect.</param>
        /// <param name="positionDependence">The position dependence of the defect.</param>
        public DefectParameters(double energy, double fwhm, double magnitude, string label, string positionDependence)
        {
            // Copy the parameters to their corresponding properties.
            Energy = energy;
            FWHM = fwhm;
            Magnitude = magnitude;
            Label = label;
            PositionDependence = positionDependence;
        }

        #endregion [Constructor]
    }
}
