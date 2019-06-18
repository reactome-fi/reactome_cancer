/*
 * Created on May 1, 2014
 *
 */
package org.reactome.cancer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.inference.TestUtils;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.MathUtilities;

/**
 * @author gwu
 *
 */
public class PancreaticCancerAnalyzerForXu {
    
    /**
     * Default constuctor.
     */
    public PancreaticCancerAnalyzerForXu() {
    }
    
    @Test
    public void extractLFNGExpressionValues() throws Exception {
        String dir = "/Users/gwu/datasets/TCGA/PAAD/gdac.broadinstitute.org_PAAD.mRNAseq_Preprocess.Level_3.2014041600.0.0/";
        String file = dir + "PAAD.uncv2.mRNAseq_RSEM_normalized_log2.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(file);
        String line = fu.readLine();
        String[] headers = line.split("\t");
        String[] genes = new String[] {
                "LFNG",
                "RB1",
                "NOTCH3",
                "HES1",
                "HEY1"
        };
        List<String> geneList = Arrays.asList(genes);
        Map<String, List<Double>> geneToValues = new HashMap<String, List<Double>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gene = tokens[0].split("\\|")[0];
            if (geneList.contains(gene)) {
                List<Double> values = new ArrayList<Double>();
                for (int i = 1; i < tokens.length; i++)
                    values.add(new Double(tokens[i]));
                geneToValues.put(gene, values);
            }
        }
        fu.close();
        List<Double> lfngValues = geneToValues.get("LFNG");
        System.out.println("Sample\tLFNG_Exp");
        Map<String, Double> sampleToValue = new HashMap<String, Double>();
        List<String> normalSamples = new ArrayList<String>();
        for (int i = 1; i < headers.length; i++) {
            sampleToValue.put(headers[i], lfngValues.get(i - 1));
            if (headers[i].endsWith("-11")) {
                // Normal samples
                normalSamples.add(headers[i]);
                continue;
            }
            String patient = headers[i];
            int index = patient.lastIndexOf("-");
            patient = patient.substring(0, index);
            System.out.println(patient + "\t" + lfngValues.get(i - 1));
        }
        // Do a simple t-test
        System.out.println("\nTotal normal samples: " + normalSamples.size());
        double[] normalValues = new double[normalSamples.size()];
        double[] cancerValues = new double[normalSamples.size()];
        for (int i = 0; i < normalSamples.size(); i++) {
            String normalSample = normalSamples.get(i);
            String cancerSample = normalSample;
            int index = cancerSample.lastIndexOf("-");
            cancerSample = cancerSample.substring(0, index) + "-01";
            normalValues[i] = sampleToValue.get(normalSample);
            cancerValues[i] = sampleToValue.get(cancerSample);
        }
        
        System.out.println("Normal values: " + outputArray(normalValues));
        System.out.println("Cancer values: " + outputArray(cancerValues));
        System.out.println("Pair-wise t: " + TestUtils.pairedT(normalValues, cancerValues));
        System.out.println("Pair-wise t-test: " + TestUtils.pairedTTest(normalValues, cancerValues));
    
        System.out.println("\nCorrelation containing 3 normal samples:");
        List<Double> rb1Values;
        calculateCorrelations(geneToValues);
        System.out.println("\nCorrelation after removing 3 normal samples:");
        for (String gene : geneList) {
            List<Double> values = geneToValues.get(gene);
            filterOutNormalValues(headers, values);
        }
        calculateCorrelations(geneToValues);
    }

    private void calculateCorrelations(Map<String, List<Double>> geneToValues) throws MathException {
        // Calculate correlation
        List<Double> rb1Values = geneToValues.get("RB1");
        for (String gene : geneToValues.keySet()) {
            List<Double> values = geneToValues.get(gene);
            PearsonsCorrelation corr = MathUtilities.constructPearsonCorrelation(rb1Values, values);
            System.out.println("Correlation between RB1 and " + gene + ": " +
                                corr.getCorrelationMatrix().getEntry(0, 1) + " (" +
                                corr.getCorrelationPValues().getEntry(0, 1) + ")");
        }
    }
    
    private void filterOutNormalValues(String[] headers, List<Double> values) {
        int i = 0;
        for (Iterator<Double> it = values.iterator(); it.hasNext();) {
            it.next();
            if (headers[i].endsWith("-11"))
                it.remove();
            i++;
        }
    }
    
    private String outputArray(double[] array) {
        StringBuilder builder = new StringBuilder();
        for (double value : array)
            builder.append(value).append(",");
        builder.delete(builder.length() - 1, builder.length());
        return builder.toString();
    }
    
}
