using System;
using System.IO;
using System.Linq;
using System.Reflection;
using Thermo.Magellan.BL.Data;
using Thermo.Magellan.BL.Processing;
using Thermo.Magellan.BL.Processing.Interfaces;
using Thermo.Magellan.MassSpec;
using Thermo.Magellan.Utilities;
using log4net;
using System.Collections.Generic;

namespace DiscovererNodes
{
    /// <summary>
    /// ComplementaryFinder is an implementation of Complementary Finder algorithm.
    /// <para>
    /// Note: Please, reference to Kryuchkov et al J. Proteome Res., 2013, 12(7), 3362–3371; DOI: 10.1021/pr400210m for details
    /// </para>
    /// </summary>
    [ProcessingNode("4799A783-EE7C-4F90-92C3-E495AB17AE05",
        Category = ProcessingNodeCategories.DataProcessing,
        DisplayName = "Complementary Finder",
        Description = @"Identifies complementary peaks in the spectrum, performs extraction of co-eluted peptide spectra,
based on complementary pairs, performs spectra intensification. Details of the algorithm can be found in
Kryuchkov et al. J. Proteome Res., 2013, 12(7), 3362–3371; DOI: 10.1021/pr400210m",
        MainVersion = 1,
        MinorVersion = 1)]
    [ProcessingNodeConstraints(UsageConstraint = UsageConstraint.Unrestricted)]
    [ConnectionPoint(
        "Incoming",
        ConnectionDirection = ConnectionDirection.Incoming,
        ConnectionMultiplicity = ConnectionMultiplicity.Single,
        ConnectionMode = ConnectionMode.Manual,
        ConnectionRequirement = ConnectionRequirement.RequiredAtDesignTime,
        ConnectionDataHandlingType = ConnectionDataHandlingType.InMemory)]
    [ConnectionPointDataContract(
        "Incoming",
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
        MassSpecDataTypes.MSnSpectra)]

    public class ComplementaryFinder : ProcessingNode<MassSpectrumCollection>, IResultsSink<MassSpectrumCollection>
        /* Encapsulate Complementary Finder funcionality */
    {
        private const int PackageSize = 1000;
        private MassSpectrumCollection spectrumBuffer = new MassSpectrumCollection(PackageSize); //send buffer
        private static readonly ILog Log = LogManager.GetLogger(MethodBase.GetCurrentMethod().DeclaringType);

        const double Proton = 1.00727646677;

        private delegate MassCentroid PeakIntensifier(MassCentroid peak, double rank); 
        private PeakIntensifier intensify = null;
        private MultiRange exclusionList = null;

        private int spectra_in = 0;
        private int spectra_out = 0;
        private int magic = 0;

        [IntSelectionParameter(
            Category = "1. Spectrum extraction",
            DisplayName = "Secondary peptide charge",
            Description = "Possible charge states of secondary peptides",
            SelectionValues = new int[] { 1, 2, 3, 4, 5, 6, 7, 8 },
            DefaultValue = "2, 3, 4",
            IsMultiSelect = true,
            ValueRequired = true,
            Position = 1)]

        public SimpleSelectionParameter<int> chargeStates;

        [IntegerParameter(
            Category = "1. Spectrum extraction",
            DisplayName = "Min number of peaks",
            Description = "Minimum number of peaks per secondary spectrum to extract it",
            MinimumValue = "2",
            DefaultValue = "6",
            ValueRequired = true,
            Position = 2)]

        public IntegerParameter minPeaks;

        [MassToleranceParameter(
            Category = "1. Spectrum extraction",
            DisplayName = "Parent mass tolerance",
            Description = "Mass tolerance to consider parent masses the same",
            MinimumValue = "0.0 ppm",
            MaximumValue = "100.0 ppm",
            DefaultValue = "5.0 ppm",
            Subset = "ppm",
            ValueRequired = true,
            Position = 3)]

        public MassToleranceParameter mergePPM;

        [BooleanParameter(
            Category = "1. Spectrum extraction",
            DisplayName = "Relative borders",
            Description = "If True the lower and upper borders of allowed coisolation window (below) are set relative to the instrument isolation window, absolute value is used otherwise",
            ValueRequired = true,
            DefaultValue = "True",
            IsAdvanced = false,
            Position = 4)]

        public BooleanParameter relativeBorders;

        [DoubleParameter(
            Category = "1. Spectrum extraction",
            DisplayName = "Lower border of coisolation window (u)",
            Description = @"Lower border of the mass window to consider coisolated peptides.
Positive values move the border to lower mass values (increase isolation window). Negative values should be used only with relative borders.
If not set the isolation window saved in MS spectrum is used",
            MinimumValue = "-10.0",
            MaximumValue = "10.0",
            DefaultValue = "0.6",
            ValueRequired = false,
            IsAdvanced = false,
            Position = 5)]

        public DoubleParameter lowerBorder;

        [DoubleParameter(
            Category = "1. Spectrum extraction",
            DisplayName = "Upper border of coisolation window (u)",
            Description = @"Upper border of the mass window to consider coisolated peptides.
Positive values move the border to higher mass values (increase isolation window). Negative values should be used only with relative borders.
If not set the isolation window saved in MS spectrum is used",
            MinimumValue = "-10.0",
            MaximumValue = "10.0",
            ValueRequired = false,
            IsAdvanced = false,
            Position = 6)]

        public DoubleParameter upperBorder;

        [StringSelectionParameter(
            Category = "1. Spectrum extraction",
            DisplayName = "Spectra to return",
            Description = "Which spectra should be returned for further processing",
            SelectionValues = new string[] { "Primary and Secondary", "Primary Only", "Secondary Only", "Intensify Only" },
            DefaultValue = "Primary and Secondary",
            IsMultiSelect = false,
            ValueRequired = true,
            IsAdvanced = true,
            Position = 7)]

        public SimpleSelectionParameter<string> returnSelection;

        [BooleanParameter(
            Category = "4. Secondary spectra mass verification",
            DisplayName = "Enable MS1 verification",
            Description = "If enabled checks if the predicted mass peak for secondary spectra is present in parent MS1 spectrum and removes spectra that don't meet this criteria",
            DefaultValue = "False",
            ValueRequired = true,
            IsAdvanced = true,
            Position = 1)]

        public BooleanParameter verifyMS1;

        [MassToleranceParameter(
            Category = "4. Secondary spectra mass verification",
            DisplayName = "Mass tolerance",
            Description = "Mass tolerance used for MS1 verification",
            MinimumValue = "0.0 mmu| 0.0 ppm",
            MaximumValue = "100.0 mmu| 100.0 ppm",
            DefaultValue = "5.0 ppm",
            Subset = "mmu|ppm",
            ValueRequired = true,
            IsAdvanced = true,
            Position = 2)]

        public MassToleranceParameter VerifyPPM;

        [BooleanParameter(
            Category = "2. Peak intensification",
            DisplayName = "Enable intensification",
            Description = "Set true if peak intensification is necessary",
            DefaultValue = "True",
            ValueRequired = true,
            Position = 1)]

        public BooleanParameter doIntensification;

        [IntegerParameter(
            Category = "2. Peak intensification",
            DisplayName = "Intensification factor: b - y",
            Description = "Intensifiaction factor for b - y ion pairs",
            DefaultValue = "1",
            ValueRequired = true,
            Position = 2)]

        public IntegerParameter byRank;

        [StringParameter(
            Category = "3. Peak exclusion",
            DisplayName = "Exclusion list",
            Description = "Path to the exclusion list (the default one contains all immonium ions > 50 Da)",
            DefaultValue = ".\\cf.exclude",
            ValueRequired = true,
            Position = 1)]

        public StringParameter ignoreListPath;

        [MassToleranceParameter(
            Category = "3. Peak exclusion",
            DisplayName = "Exclusion tolerance",
            Description = "Exclusion window is the doubled tolerance, e.g. (-5 .. +5) ppm",
            MinimumValue = "0.0 mmu| 0.0 ppm",
            MaximumValue = "100.0 mmu| 100.0 ppm",
            DefaultValue = "10.0 ppm",
            Subset = "mmu|ppm",
            ValueRequired = true,
            Position = 2)]

        public MassToleranceParameter ignoreListPPM;

        private MassCentroid validIntensify(MassCentroid peak, double maxIntensity)
        /* Change intensity of peak according to rank */
        {
            MassCentroid newPeak = peak.Clone();

            newPeak.Intensity += this.byRank.Value * maxIntensity;

            return newPeak;
        }

        private MassCentroid voidIntensify(MassCentroid peak, double maxIntensity)
        /* Void intensification */
        {
            return peak;
        }

        private short getCharge (int index)
        /* Translate index of MultiRange to charge state */
        {
            return (short) this.chargeStates.Values[index / 2];
        }

        private MultiRange loadExcludes ()
        {
            MultiRange exclusionList = new MultiRange();

            using (StreamReader massinput = new StreamReader(ignoreListPath.Value))
            {
                double mass, ppm;
                do
                {
                    mass = Double.Parse(massinput.ReadLine());
                    ppm = ignoreListPPM.Value.GetTolerance(mass, MassToleranceUnit.U);
                    exclusionList.Add(mass - ppm, mass  + ppm);
                } while (!massinput.EndOfStream);
            }

            return exclusionList;
        }
        
        private MassCentroidCollection annotateSpectrum(MassSpectrum spectrum, MultiRange range)
        /* Search for peak pairs that reside inside MultiRange
         * Annotate b - y pairs in the spectrum with unique parent mass ID
         * Return: MassCentroidCollection of all parent masses
         */
        {
            int start = 0; //lower seeker
            int end = spectrum.PeakCentroids.Count - 1; //upper seeker
            int k; //intermediate seeker
            int lowIndex; //position in multirange
            double mass;
            MassAggregator masses = new MassAggregator(this.mergePPM.Value.Tolerance); //collection of possible masses
            int spectrumID; //unique ID for each subspectrum (set of ion pairs) in spectrum

            MassCentroid instrumentMass = spectrum.Precursor.MeasuredMonoisotopicPeakCentroids[0].Clone();
            instrumentMass.Charge = spectrum.Precursor.Charge;
            instrumentMass.Position = (instrumentMass.Position - Proton) * instrumentMass.Charge;
            instrumentMass.SignalToNoise = 0;

            masses.AddCentroid(instrumentMass); //always add initial mass

            while (start < end) //continue untill two seekers are met
            {
                mass = spectrum.PeakCentroids[start].Position + spectrum.PeakCentroids[end].Position - 2 * Proton; //sum of two peaks converted to neutral M

                if (mass > range.Highest) //mass is too big
                {
                    end--; //make smaller = move end backwards
                }
                else if (mass < range.Lowest) //mass is too small
                {
                    start++; //make bigger = move start forward
                }
                else
                {
                    //check all possible pairs with masses higer or equal than current 
                    k = 0;
                    lowIndex = 0;

                    do
                    {
                        lowIndex = range.lastIndexBelow(lowIndex, mass); //find the highest index in range smaller than mass
                        if (lowIndex % 2 == 0) //if *lowIndex* is even (it is the begginning of subrange) add peak
                        {
                            spectrumID = masses.AddCentroid(
                                new MassCentroid(mass, spectrum.PeakCentroids[start + k].Intensity + spectrum.PeakCentroids[end].Intensity,
                                    getCharge(lowIndex), spectrum.Precursor.Resolution, 2));
                            //signal-to-noise is the number of peaks supporting this mass
                            spectrum.PeakCentroids[start + k].SetAnnotation(spectrumID.ToString(), true);
                            spectrum.PeakCentroids[end].SetAnnotation(spectrumID.ToString(), true);
                            k++; //try next peak
                        }
                        else //find the index of the first peak above the beginning of next subrange
                        {
                            k = spectrum.PeakCentroids.FindIndex(start + k, //start from current peak
                                x => x.Position > spectrum.PeakCentroids[end].Position //peaks above [end] were checked earlier
                                    || x.Position > range[lowIndex + 1] - spectrum.PeakCentroids[end].Position + 2 * Proton) - start;
                        }

                        //Console.WriteLine("start: {0}, k: {1}, end: {2}", start, k, end);
                        if (start + k >= end || k < 0) break; // no such peak

                        mass = spectrum.PeakCentroids[start + k].Position + spectrum.PeakCentroids[end].Position - 2 * Proton;

                    } while (mass < range.Highest && k < spectrum.PeakCentroids.Count - 1);
                    //do so untill mass exeeds the upper border or there is no peaks in the spectrum

                    end--; // make initial mass smaller, hiher masses were checked earlier
                }
            }

            //Log.DebugFormat("Annotation: Scan {0}; {1} peaks; {2} possible masses", spectrum.Header.ScanNumbers[0], spectrum.PeakCentroids.Count, masses.Count);

            return masses.getAll();
        }

        private MassSpectrumCollection extractSpectra(MassSpectrum spectrum, MassCentroidCollection parentMasses, int minPeaks)
            /* Extract spectra according to annotation
             * Only parent masses supported by at least minPeaks peaks are used,
             * parentMasses is a collection of masses
             * Return: MassSpectrumCollection of spectra, position 0 - primary, other - secondary
             */
        {
            MassSpectrumCollection outSpectra = new MassSpectrumCollection();
            List<int> validIDs = new List<int>();

            validIDs.Add(0); //always include the initial mass (index 0)

            if (verifyMS1.Value) //verify parent masses according to MS1 if necessary
            {
                foreach (int i in areMassesInSpectrum(parentMasses, spectrum)) //check IDs, supported by MS1
                {
                    if (parentMasses[i].SignalToNoise >= minPeaks) validIDs.Add(i); //signal-to-noise is used to store number of ions in a group
                }
            }
            else
            {
                for (int i = 1; i < parentMasses.Count; i++) //check all IDs 
                {
                    if (parentMasses[i].SignalToNoise >= minPeaks) validIDs.Add(i);
                }
            }

            outSpectra.Add(spectrum.Clone()); //clone initial spectrum to position 0; all unassigned peaks will be here

            //Log.DebugFormat("Extraction: Scan {0}; {1} valid mass(es)", spectrum.Header.ScanNumbers[0], validIDs.Count);

            if (validIDs.Count > 0) //if at least one valid spectrum was found
            {
                for (int i = 0; i < validIDs.Count; i++)//add one clone for each additional spectrum
                {
                    MassSpectrum additionalSpectrum = spectrum.Clone();

                    additionalSpectrum.ProfilePoints.Clear(); //remove all peaks
                    additionalSpectrum.PeakCentroids.Clear();
                    additionalSpectrum.Precursor.Charge = parentMasses[validIDs[i]].Charge; //correct mass
                    additionalSpectrum.Precursor.MeasuredMonoisotopicPeakCentroids[0] = parentMasses[validIDs[i]];
                    additionalSpectrum.Precursor.SinglyChargedMass = new MassValue(parentMasses[validIDs[i]].Position + Proton);
                    additionalSpectrum.Precursor.MeasuredMonoisotopicPeakCentroids[0].Position = (parentMasses[validIDs[i]].Position //calculate m/z
                        + parentMasses[validIDs[i]].Charge * Proton) / parentMasses[validIDs[i]].Charge;
                    //?correct scan number?
                    outSpectra.Add(additionalSpectrum);
                }

                foreach (MassCentroid peak in spectrum.PeakCentroids) //try each Centroid in initial spectrum
                {
                    if (peak.HasAnnotations) //if it was annotated at all (most of peaks aren't annotated)
                    {
                        for (int i = 0; i < validIDs.Count; i++) //check if it was annotated with valid ID
                        {
                            if (peak.GetAnnotation(validIDs[i].ToString()) != null)
                            {
                                outSpectra[i + 1].PeakCentroids.Add(intensify(peak, spectrum.Header.BasePeakIntensity)); //add it to the corresponding outSpectra element
                                outSpectra[0].PeakCentroids.Remove(peak); //remove from initial
                            }
                        }
                    }
                }

                //Log.DebugFormat("Extraction: Subspectrum 0 - {0} peaks", outSpectra[0].PeakCentroids.Count);

                for (int i = 1; i < outSpectra.Count; i++) //add unassigned peaks to each spectrum
                {
                    //Log.DebugFormat("Extraction: Subspectrum {0} - {1} peaks", i, outSpectra[i].PeakCentroids.Count);
                    outSpectra[i].PeakCentroids.AddRange(outSpectra[0].PeakCentroids.AsEnumerable<MassCentroid>());
                    outSpectra[i].PeakCentroids.Sort();
                }

                outSpectra.RemoveAt(0); //remove the initial spectrum
            }

            return outSpectra;
        }

        private MassSpectrumCollection intensifyOnly(MassSpectrum spectrum)
            /* Do intensification of peaks corresponding to targeted mass only, all secondary masses are ignored
             * spectrum - annotated spectrum
             * Return: MassSpectrumCollection with one spectrum
             */
        {
            MassSpectrumCollection outSpectra = new MassSpectrumCollection(1);
            MassSpectrum unintSpectrum = spectrum.Clone(); //position 0 is for unintensified peaks

            MassSpectrum intSpectrum = spectrum.Clone(); //create new entity of spectrum
            intSpectrum.ProfilePoints.Clear(); //remove all peaks
            intSpectrum.PeakCentroids.Clear();

            foreach (MassCentroid peak in spectrum.PeakCentroids) //try each Centroid in initial spectrum
            {
                if (peak.HasAnnotations) //if it was annotated at all (most of peaks aren't annotated)
                {
                    if (peak.GetAnnotation("0") != null) //if the peak is annotated to the targeted mass
                    {
                        intSpectrum.PeakCentroids.Add(intensify(peak, spectrum.Header.BasePeakIntensity)); //add intensified peak
                        unintSpectrum.PeakCentroids.Remove(peak); //remove from unintensified peaks
                    }
                }
            }

            intSpectrum.PeakCentroids.AddRange(unintSpectrum.PeakCentroids.AsEnumerable<MassCentroid>()); //add all unintensified peaks
            intSpectrum.PeakCentroids.Sort();

            outSpectra.Add(intSpectrum);//add spectrum to the output
            return outSpectra;
        }

        private MassCentroidCollection applyExclusionList (MassSpectrum spectrum)
            /* Process the spectrum (IT WILL BE MODIFIED!) by removing all peaks that are inside exclusion list
             * collect all filtered peaks in a separate collection and return it
             */
        {
            MassCentroidCollection filtered = new MassCentroidCollection(); //filtered peaks collection

            for (int i = 0; i < exclusionList.Count; i += 2) //check all ranges in exclusion list
            {
                foreach (MassCentroid peak in spectrum.PeakCentroids.FindAllPeaksWithinMassRange(exclusionList[i], exclusionList[i + 1]))
                    //find all peaks inside the range; remove them from the spectrum and save in collection
                {
                    spectrum.PeakCentroids.Remove(peak);
                    filtered.Add(peak);
                }
            }

                return filtered;
        }
        
        private MassSpectrumCollection rejoinSpectra (MassSpectrumCollection inSpectra, MassCentroidCollection peaksToJoin)
            /* Add peaks from MassCentroidCollection to every spectrum in MassSpectrumCollection
             * return modified MassSpectrumCollection
             */
        {
            foreach (MassSpectrum spectrum in inSpectra)
            {
                spectrum.PeakCentroids.AddRange(peaksToJoin.AsEnumerable());
                spectrum.PeakCentroids.Sort();
            }

            return inSpectra;
        }

        private List<int> areMassesInSpectrum(MassCentroidCollection masses, MassSpectrum spectrum)
            /* For each mass in mass list search if the peak is present in the spectrum
             * Collect and return indices of valid peaks */
        {
            double tol; //tolerance in amu
            double mz, cmz; //predicted m/z and the closest match from MS1
            List<int> valid = new List<int>();

            for (int c = 1; c < masses.Count; c++) //check each mass except the first one (first mass - primary spectrum)
            {
                mz = (masses[c].Position + Proton * masses[c].Charge) / masses[c].Charge; //calculate m/z value to search and mass tolerance
                tol = this.VerifyPPM.Value.GetToleranceInU(mz);

                //check if there is at least one peak in precursor spectrum and find the closest peak to m/z value
                if (spectrum.Precursor.IsotopeClusterPeakCentroids.Count > 0) cmz = spectrum.Precursor.IsotopeClusterPeakCentroids.FindClosestPeak(mz).Position;
                else cmz = 0;

                if (Math.Abs(cmz - mz) <= tol) valid.Add(c); //add valid index
            }

            return valid;
        }

        private void ProcessBatch(MassSpectrumCollection results)
        /* Process a set of mass spectra using complementary finder */
        {
            foreach (MassSpectrum iSpec in results)
            {
                double low, high;
                MultiRange range = new MultiRange();
                MassSpectrum spectrum = iSpec.Clone(); //make local copy of the input spectrum

                if (!this.lowerBorder.IsValueSet) low = spectrum.ScanEvent.IsolationMass - spectrum.ScanEvent.IsolationWindow.LowerLimit;//get right value of lower border
                else if(this.relativeBorders.Value) low = spectrum.ScanEvent.IsolationMass - spectrum.ScanEvent.IsolationWindow.LowerLimit + this.lowerBorder.Value;
                else low =  this.lowerBorder.Value;

                if (!this.upperBorder.IsValueSet) high = spectrum.ScanEvent.IsolationWindow.UpperLimit - spectrum.ScanEvent.IsolationMass;//get right value of upper border
                else if (this.relativeBorders.Value) high = spectrum.ScanEvent.IsolationWindow.UpperLimit - spectrum.ScanEvent.IsolationMass + this.upperBorder.Value;
                else high =  this.upperBorder.Value;

                //mgf-file input compatability
                if (spectrum.ScanEvent.IsolationMass == 0.0) spectrum.ScanEvent.IsolationMass = spectrum.Precursor.MeasuredMonoisotopicPeakCentroids[0].Position;
                spectrum.PeakCentroids.Sort();

                foreach (int z in this.chargeStates.Values)
                {
                    range.Add((spectrum.ScanEvent.IsolationMass - low - Proton) * z, (spectrum.ScanEvent.IsolationMass + high - Proton) * z); //neutral MolMass
                }

                if (range.Count == 0)//check if parameters are valid
                {
                    SendAndLogErrorMessage("The window for possible parent masses is empty! Please, check the settings of allowed coisolation window", true);
                    throw new Exception("Allowed coisolation window is empty");
                }

                //for debugging purposes
                if (spectrum.Header.ScanNumbers[0] == magic)
                {
                    Log.DebugFormat("Magic!");
                }

                MassCentroidCollection excludedPeaks = applyExclusionList(spectrum); //apply exclusion list and put all excluded peaks in separate collection

                MassCentroidCollection parentMasses = annotateSpectrum(spectrum, range); //annotate complementary pairs

                MassSpectrumCollection extractedSpectra = new MassSpectrumCollection();

                switch (this.returnSelection.Value) //filter spectra acording to return selection
                {
                    case "Primary and Secondary":
                        extractedSpectra = rejoinSpectra(extractSpectra(spectrum, parentMasses, this.minPeaks.Value), excludedPeaks); //extract spectra
                        break; 
                    case "Primary Only":
                        extractedSpectra = rejoinSpectra(extractSpectra(spectrum, parentMasses, this.minPeaks.Value), excludedPeaks); //extract spectra
                        extractedSpectra.RemoveRange(1, extractedSpectra.Count - 1);
                        break; //remove all secondary
                    case "Secondary Only":
                        extractedSpectra = rejoinSpectra(extractSpectra(spectrum, parentMasses, this.minPeaks.Value), excludedPeaks); //extract spectra
                        extractedSpectra.RemoveAt(0);
                        break; //remove primary
                    case "Intensify Only":
                        extractedSpectra = rejoinSpectra(intensifyOnly(spectrum), excludedPeaks); 
                        break;//just intensify ions related to targeted mass
                    default:
                        throw new Exception("Unrecognized return selection parameter"); //control
                }

                toSend(extractedSpectra); //transfer spectra to sending buffer
            }
        }

        public void OnResultsSent(IProcessingNode sender, MassSpectrumCollection result)
        {
            if (intensify == null) //initialize intensification method
            {
                if (this.doIntensification.Value) intensify = validIntensify;
                else intensify = voidIntensify;
            }

            if (exclusionList == null) //read exclusion list
            {
                try
                {
                    exclusionList = loadExcludes();
                }
                catch (Exception ex)
                {
                    SendAndLogErrorMessage(String.Format("Can't load exlusion list from '{0}': {1}", ignoreListPath.Value, ex.Message));
                    throw;
                }
            }

            Log.Debug(String.Format("{0} spectra received", result.Count));
            spectra_in += result.Count;

            try
            {
                ProcessBatch(result);
            }
            catch (Exception ex)
            {
                SendAndLogErrorMessage("Error: " + ex.Message);
                throw;
            }
        }

        public override void OnParentNodeFinished(IProcessingNode sender, ResultsArguments eventArgs)
        {
            flushSend();
            FireProcessingFinishedEvent(new ResultsArguments(MassSpecDataTypes.MSnSpectra));
        }

        private void toSend (MassSpectrumCollection inspectra)
        /* Check if results can still fit in buffer and add them, otherwise send buffer and flush it */
        {
            //Log.DebugFormat("toSend received {0} spectra", inspectra.Count);

            ProcessingServices.SpectrumProcessingService.InitializeSpectra(this, inspectra);
            if (inspectra.Count + spectrumBuffer.Count > PackageSize)
            {
                SendResults(spectrumBuffer);

                Log.DebugFormat(String.Format("{0} spectra sent", spectrumBuffer.Count));
                SendAndLogTemporaryMessage(String.Format("{0} spectra received, {1} spectra sent", spectra_in, spectra_out), true);
                spectra_out += spectrumBuffer.Count;

                spectrumBuffer = new MassSpectrumCollection(PackageSize);
            }

            spectrumBuffer.AddRange(inspectra.AsEnumerable());
        }

        private void flushSend()
        /* Flush the all spectra in buffer */
        {
            SendResults(spectrumBuffer);
            Log.DebugFormat(String.Format("{0} spectra sent", spectrumBuffer.Count));
            spectra_out += spectrumBuffer.Count;
            SendAndLogMessage(String.Format("Processing finished: {0} spectra received, {1} spectra sent", spectra_in, spectra_out), true);
        }
    }

    class MultiRange : Object
    /* Represents an range consisting of several subranges
     * eg. (1.0, 2.1) and (5.02, 50.04) and ...
     * NOTE: borders aren't included
     */
    {
        private List<double> points = new List<double>();

        public double Lowest //lowest border of all ranges
        {
            get
            {
                if (points.Count > 0) return points[0];
                else throw new Exception("Empty MultiRange");
            }
        }

        public double Highest //highest border of all ranges
        {
            get
            {
                if (points.Count > 0) return points[points.Count - 1];
                else throw new Exception("Empty MultiRange");
            }
        }

        public MultiRange()
        {
            //empty
        }

        public void Add(double lower, double upper)
        //add a range to multirange
        {
            if (lower >= upper) return; //length of the range should be positive

            int first = 0; //position to insert values in points list
            int last = 0;

            while (first < points.Count && lower > points[first]) //find position of lower
            {
                first++;
            }

            last = first;

            while (last < points.Count && upper > points[last]) //find position of upper
            {
                last++;
            }

            if (first % 2 == 0 && last % 2 == 0) // complete pair to the left and to the right
            {
                points.RemoveRange(first, last - first);
                points.Insert(first, upper);
                points.Insert(first, lower);
            }
            else if (first % 2 == 0 && last % 2 == 1) //complete pair to the left, incomplete to the right
            {
                points.RemoveRange(first, last - first);
                points.Insert(first, lower);
            }
            else if (first % 2 == 1 && last % 2 == 0) //incomplete pair to the left, complete to the right
            {
                points.RemoveRange(first, last - first);
                points.Insert(first, upper);
            }
            else //incomplete pair to the left and to the right
            {
                points.RemoveRange(first, last - first);
            }

        }

        public bool isIn(double value)
        //check if value is in multirnage
        {
            if (points.Count % 2 == 1) throw new Exception("Invalid range"); //internal check

            for (int i = 0; i < points.Count; i += 2)
            {
                if (value > points[i] && value < points[i + 1]) return true;
            }

            return false;
        }

        public int lastIndexBelow(int start, double value)
        {
            return points.FindIndex(start, x => x > value) - 1;
        }

        public int firstIndexAbove(int start, double value)
        {
            return points.FindIndex(start, x => x > value);
        }

        public double this[int index]
        {
            get { return points[index]; }
        }

        public int Count
        {
            get 
            {
                return points.Count;
            }
        }
    }

    class MassAggregator : Object
    /* Encloses list of MassCentroids, that are kept unique and each has a unique ID
     * Centroids withing certian precision are merged on addition and total intensity is stored
     */
    {
        private MassCentroidCollection massList = new MassCentroidCollection();
        private double ppm;

        public MassAggregator(double precision)
        /* Initialize with certain ppm precision */
        {
            ppm = precision * 1e-6;
        }

        public int AddCentroid(MassCentroid newMass)
        /* Check if mass can be merged with some of existing and if so the ID is returned; 
         * otherwise mass is added to the list and new ID is generated */
        {
            for (int i = 0; i < massList.Count; i++) //check all masses, that are already in the list
            {
                if (Math.Abs(massList[i].Position - newMass.Position) < massList[i].Position * ppm)
                {
                    //merge peaks
                    massList[i].Position = (massList[i].Position * massList[i].Intensity + newMass.Position * newMass.Intensity) /
                        (massList[i].Intensity + newMass.Intensity);
                    massList[i].Intensity = massList[i].Intensity + newMass.Intensity;
                    massList[i].SignalToNoise = massList[i].SignalToNoise + newMass.SignalToNoise;

                    return i;
                }
            }
            massList.Add(newMass);
            return massList.Count - 1;
        }

        public MassCentroid this[int index]
        {
            get { return massList[index]; }
        }

        public int Count
        {
            get { return massList.Count; }
        }

        public MassCentroidCollection getAll()
        {
            return massList;
        }
    }

    
}
