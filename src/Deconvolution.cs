using System;
using System.Collections.Generic;
using System.Reflection;
using System.Linq;
using Thermo.Magellan.BL.Data;
using Thermo.Magellan.BL.Processing;
using Thermo.Magellan.BL.Processing.Interfaces;
using Thermo.Magellan.MassSpec;
using Thermo.Magellan.Utilities;
using System.Xml.Serialization;
using log4net;

namespace DiscovererNodes
{
    /// <summary>
    /// Deconvolution node performs MSn spectrum deconvolution.
    /// </summary>
    [ProcessingNode("30C83F6B-F63C-49D7-B6EE-17671BBC5A99",
        Category = ProcessingNodeCategories.SpectrumProcessing,
        DisplayName = "MSn Deconvolution",
        Description = "Performs deconvolution/decharging/deisotoping of spectra.",
        MainVersion = 1,
        MinorVersion = 1)]
    [ProcessingNodeConstraints(UsageConstraint = UsageConstraint.Unrestricted)]
    [ConnectionPoint(
        "IncomingSpectra",
        ConnectionDirection = ConnectionDirection.Incoming,
        ConnectionMultiplicity = ConnectionMultiplicity.Single,
        ConnectionMode = ConnectionMode.Manual,
        ConnectionRequirement = ConnectionRequirement.RequiredAtDesignTime,
        ConnectionDisplayName = ProcessingNodeCategories.SpectrumAndFeatureRetrieval,
        ConnectionDataHandlingType = ConnectionDataHandlingType.InMemory)]
    [ConnectionPointDataContract(
        "IncomingSpectra",
        MassSpecDataTypes.MSnSpectra)]
    [ConnectionPoint(
        "OutgoingSpectra",
        ConnectionDirection = ConnectionDirection.Outgoing,
        ConnectionMultiplicity = ConnectionMultiplicity.Multiple,
        ConnectionMode = ConnectionMode.Manual,
        ConnectionRequirement = ConnectionRequirement.Optional,
        ConnectionDataHandlingType = ConnectionDataHandlingType.InMemory)]
    [ConnectionPointDataContract(
        "OutgoingSpectra",
        MassSpecDataTypes.MSnSpectra,
        DataTypeAttributes = new[] { MassSpecDataTypeAttributes.Decharged, MassSpecDataTypeAttributes.Deconvoluted, MassSpecDataTypeAttributes.Deisotoped })]
    public class Deconvolution : SpectrumProcessingNode
    {
        const double Proton = 1.00727646677; // Mass of proton
        private static readonly ILog Log = LogManager.GetLogger(MethodBase.GetCurrentMethod().DeclaringType);

        IsotopePattern findGroup(MassSpectrum spectrum, MassCentroid peak, ref int id)
        {
            double isotopediff = IsoDiff.Value.Mass;
            double isotopeppm = IsoPPM.Value.Tolerance;

            var peaks = new MassCentroidCollection(); // Centroids of isotopic cluster
            double last = peak.Position; // position of the last peak in the cluster

            peaks.Add(peak); // add the starting peak anyway and remove it from spectrum
            spectrum.PeakCentroids.Remove(peak); 

            while (true) // try untill the end of spectrum
            {
                var nearest = spectrum.PeakCentroids.FindClosestPeak(last + isotopediff / peak.Charge);//find the peak nearest to the predicted one

                if (nearest != null && (Math.Abs(nearest.Position - last - isotopediff / peak.Charge) < nearest.Position * isotopeppm / 1000000))
                // null means - there are no peaks anymore; if the difference is under the precision border add the peak to cluster
                {
                    peaks.Add(nearest);
                    last = nearest.Position; // redefine last peak in cluster
                    spectrum.PeakCentroids.Remove(nearest); // remove the peak found from the spectrum
                }
                else //  return the result
                {
                    IsotopePattern group = new IsotopePattern();
                    group.PatternID = id++; // increment the ID
                    group.PatternPeaks = peaks; // add peaks, resolution and charge
                    group.Charge = peak.Charge;
                    group.Resolution = peak.Resolution;

                    return group;
                }
            }
        }

        void mergePeaks(MassCentroidCollection spectrum)
        {
            //Process spectrum !in place! merging peaks that are within ppm from each other
            //helper function to process raw deconvoluted spectrum
            int i = 1; //start with the second peak and go until the end
            while (i < spectrum.Count)
            {
                if (spectrum[i].Position - spectrum[i - 1].Position < spectrum[i].Position * MergePPM.Value.Tolerance / 1e6) //peaks can be merged
                {
                    spectrum[i - 1].Position = (spectrum[i].Position + spectrum[i - 1].Position) / 2;
                    spectrum[i - 1].Intensity = spectrum[i].Intensity + spectrum[i - 1].Intensity;
                    spectrum.RemoveAt(i);
                }
                else
                    i++;
            }
        }

        List<IsotopePattern> spectrum2IPList(MassSpectrum ospectrum) // Identify isotopic patterns in the spectrum and return the list of them
        {
            List<IsotopePattern> patterns = new List<IsotopePattern>();
            IsotopePattern cluster;

            MassSpectrum spectrum = ospectrum.Clone();//clone the original spectrum, since later it will be modified
            MassCentroid peak;
            int id = 0; //ID for isotopic pattern
            int i = 0; //index

            while (i < spectrum.PeakCentroids.Count) // do untill at least one peak is in the spectrum
            {

                peak = spectrum.PeakCentroids[i]; // take the next peak

                if (peak.Charge > 0) //if the charge is assigned try finding isotopic cluster
                    patterns.Add(findGroup(spectrum, peak, ref id));
                else
                    i++; //move to next peak
            }
            if (KeepUnassigned.Value)
            {
                while (spectrum.PeakCentroids.Count > 0)//the rest should be peaks with unassigned charges
                {
                    peak = spectrum.PeakCentroids[0];
                    if (peak.Charge != 0) //just for the case
                        throw new Exception("Charge assigned peak");

                    peak.Charge = (short)DefCharge.Value; // unassigned charge goes to default
                    cluster = findGroup(spectrum, peak, ref id);
                    if (cluster.PatternPeaks.Count >= MinPeakPerClusterUA.Value) patterns.Add(cluster); // try to find isotopic peaks and add the cluster to the list
                }
            }

            return patterns;
        }

        MassCentroidCollection deconvolute(MassSpectrum spectrum) //Converts mass spectrum into sorted singly charged peak centroids
        {
            List<IsotopePattern> iplist = spectrum2IPList(spectrum);//make the list of isotopic clusters
            MassCentroidCollection spectrumPeaks = new MassCentroidCollection(); 

            foreach (IsotopePattern ip in iplist.Where(c => c.PatternPeaks.Count >= MinPeakPerClusterGlobal.Value)) // apply filter for minimal number of peaks per cluster
            {
                MassCentroid peak = new MassCentroid(ip.MonoisotopicMass * ip.Charge - (ip.Charge - 1) * Proton, ip.SummedIntensity, 1, ip.Resolution, 0); //create new centroid
                spectrumPeaks.Add(peak);
            }

            spectrumPeaks.Sort();// sort peaks
            mergePeaks(spectrumPeaks);//merge similar peaks

            return spectrumPeaks;
        }

        // Main function and parameters overriden from parent class. Functionality is here

        [BooleanParameter(
            Category = "1. Spectra deconvolution",
            DisplayName = "Keep unassigned",
            Description = "When set to True all peaks without charge assigned the default charge state",
            DefaultValue = "True",
            Position = 1)]
        public BooleanParameter KeepUnassigned;

        [IntegerParameter(
            Category = "1. Spectra deconvolution",
            DisplayName = "Default charge state",
            Description = "Default charge state for peaks without assigned charge state",
            DefaultValue = "1",
            MinimumValue = "1",
            MaximumValue = "20",
            Position = 2)]
        public IntegerParameter DefCharge;

        [MassToleranceParameter(
            Category = "1. Spectra deconvolution",
            DisplayName = "Merging peak mass tolerance",
            Description = "All peaks closer than this value after the deconvolution will be merged together.",
            DefaultValue = "5 ppm",
            MinimumValue = "0.0 ppm",
            MaximumValue = "100 ppm",
            Subset = "ppm",
            Position = 3)]
        public MassToleranceParameter MergePPM;

        [IntegerParameter(
            Category = "1. Spectra deconvolution",
            DisplayName = "Peaks per cluster (global)",
            Description = "Minimal number of isotopic peaks per cluster to process (All peaks)",
            DefaultValue = "1",
            MinimumValue = "1",
            Position = 5)]
        public IntegerParameter MinPeakPerClusterGlobal;

        [IntegerParameter(
            Category = "1. Spectra deconvolution",
            DisplayName = "Peaks per cluster (unassigned)",
            Description = "Minimal number of isotopic peaks per cluster to process (unassigned).\nWhen set to 0 it is equal to global value",
            DefaultValue = "0",
            MinimumValue = "0",
            ValueRequired = false,
            IsAdvanced = true,
            Position = 6)]
        public IntegerParameter MinPeakPerClusterUA;

        [IntegerParameter(
            Category = "1. Spectra deconvolution",
            DisplayName = "Peaks per spectrum",
            Description = "Minimal number of peaks per spectrum after the deconvolution to consider it further",
            DefaultValue = "1",
            MinimumValue = "0",
            Position = 7)]
        public IntegerParameter MinPeakPerSpectrum;

        [MassValueParameter(
           Category = "2. Isotopic cluster fitting",
           DisplayName = "Isotope peak difference",
           Description = "Mass difference between consequtive isotopic peaks\n(1.002241 u is suggested for Ga-peptides)",
           DefaultValue = "1.002713 u",
           MinimumValue = "0.1 u",
           MaximumValue = "5.0 u",
           Subset = "u",
           IsAdvanced = true,
           Position = 4)]
        public MassValueParameter IsoDiff;

        [MassToleranceParameter(
            Category = "2. Isotopic cluster fitting",
            DisplayName = "Distance tollerance",
            Description = "Mass tollerance of peak distance",
            DefaultValue = "10.0 ppm",
            MinimumValue = "0.0 ppm",
            MaximumValue = "100 ppm",
            Subset = "ppm",
            IsAdvanced = true,
            Position = 3)]
        public MassToleranceParameter IsoPPM;

        protected override MassSpectrumCollection ProcessSpectra(MassSpectrumCollection inspectra)
        {
            MassSpectrumCollection outspectra = new MassSpectrumCollection();

            if (MinPeakPerClusterUA.Value == 0) MinPeakPerClusterUA = MinPeakPerClusterGlobal;

            foreach (MassSpectrum spectrum in inspectra)
            {
                MassSpectrum dspectrum = spectrum.Clone();

                dspectrum.ProfilePoints.Clear();
                dspectrum.PeakCentroids = deconvolute(spectrum);
                
                if (dspectrum.PeakCentroids.Count >= MinPeakPerSpectrum.Value) //apply filter for minimal number of peaks per spectrum
                    outspectra.Add(dspectrum);
                
            }

            return outspectra;
        }
    }
}
