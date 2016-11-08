using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace JunctionSolver
{
    /// <summary>
    /// A class for holding physical constants.
    /// </summary>
    public static class Constants
    {
        /// <summary>
        /// The magnitude of the charge of an electron (in C).
        /// </summary>
        public static readonly double ElementaryCharge = 1.6021766e-19;

        /// <summary>
        /// The permittivity of the vacuum (in F/m).
        /// </summary>
        public static readonly double VacuumPermittivity = 8.8541878e-12;

        /// <summary>
        /// Boltzmann's constant (in J/K).
        /// </summary>
        public static readonly double BoltzmannConstant = 1.3806485e-23;
    }
}
