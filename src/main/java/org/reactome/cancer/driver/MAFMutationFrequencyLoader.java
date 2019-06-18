/*
 * Created on Apr 21, 2017
 *
 */
package org.reactome.cancer.driver;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.cancer.base.MATFileLoader;
import org.reactome.genome.Transcript;
import org.reactome.r3.UCSCDataAnalyzer;
import org.reactome.r3.util.FileUtility;

/**
 * Customized MAFFileLoader to load site mutation frequencies.
 * @author gwu
 *
 */
public class MAFMutationFrequencyLoader extends MATFileLoader {
    
    /**
     * Default constructor.
     */
    public MAFMutationFrequencyLoader() {
    }
    
    @Test
    public void testLoadGeneToNucleotideSiteMutationFrequency() throws IOException {
        String fileName = "datasets/ICGC/2016_04/Barcelona_consensus.filter.genes.dbNSFP.normalized.maf";
        Map<String, Transcript> geneToTx = new UCSCDataAnalyzer().getHumanGeneToTranscript();
        Map<String, int[]> geneToFrequencies = loadGeneToNucleotideSiteMutationFrequency(fileName, null, geneToTx);
        int[] tp53Profile = geneToFrequencies.get("TP53");
        for (int i = 0; i < tp53Profile.length; i++) {
            if (tp53Profile[i] > 0)
                System.out.println(i + "\t" + tp53Profile[i]);
        }
    }
    
    /**
     * Load the map from genes to mutation frequencies based on coding region nucleotide sequences.
     * @param fileName
     * @param excludeSamples
     * @param geneToTranscript
     * @return
     * @throws IOException
     */
    public Map<String, int[]> loadGeneToNucleotideSiteMutationFrequency(String fileName,
                                                                        Set<String> excludedSamples,
                                                                        Map<String, Transcript> geneToTranscript) throws IOException {
        Map<String, int[]> geneToSiteFrequency = new HashMap<String, int[]>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        // Headers
        String line = fu.readLine();
        String[] tokens = line.split("\t");
        int variantIndex = 0;
        int sampleIndex = 0;
        int startPositionIndex = 0;
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals("Variant_Classification"))
                variantIndex = i;
            else if (tokens[i].equals(Tumor_Sample_Barcode))
                sampleIndex = i;
            else if (tokens[i].equals("Start_position"))
                startPositionIndex = i;
        }
        Set<String> unfoundGenes = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // We will check missense_mutation only
            if (!"Missense_Mutation".equals(tokens[variantIndex]))
                continue;
            String gene = tokens[0];
            int startPosition = Integer.parseInt(tokens[startPositionIndex]);
            // It seems that there is an error in the ICGC Pancancer maf file for DCDC1
            if (gene.equals("DCDC1") && startPosition < 31284170)
                gene = "DCDC5";
            Transcript transcript = geneToTranscript.get(gene);
            if (transcript == null) {
//                throw new IllegalArgumentException("Cannot find transcript for " + gene);
//                System.err.println("Cannot find transcript for " + gene);
                unfoundGenes.add(gene);
                continue;
            }
            tokens = line.split("\t");
            String sample = parseSample(tokens[sampleIndex]);
            if (excludedSamples != null && excludedSamples.contains(sample))
                continue;
            int[] siteFrequency = geneToSiteFrequency.get(gene);
            if (siteFrequency == null) {
                Integer length = transcript.getCds().getLength();
                siteFrequency = new int[length];
                geneToSiteFrequency.put(gene, siteFrequency);
            }
            try {
                int cdsPosition = transcript.mapToCDSCoordinate(startPosition);
                if (cdsPosition <= 0) {
                    System.err.println("Mapped to negative CDS: " + cdsPosition + " in " + gene + " for " + startPosition);
                    continue;
                }
                else if (cdsPosition > siteFrequency.length - 1) {
                    System.err.println("Mapped to an outsite of CDS: " + cdsPosition + " in " + gene + " for " + startPosition);
                    continue;
                }
                siteFrequency[cdsPosition - 1] ++; // Assuming the first position is 1.
            }
            catch(IllegalArgumentException e) {
                System.err.println(e);
                continue;
            }
        }
        fu.close();
        System.err.println("Total unfound genes: " + unfoundGenes.size());
        for (String gene : unfoundGenes)
            System.err.println(gene);
        return geneToSiteFrequency;
    }
    
}
