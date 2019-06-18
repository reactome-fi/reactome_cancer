/*
 * Created on Feb 8, 2010
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.MathUtilities;

public class TCGADataHelper {

    public TCGADataHelper() {
        
    }
    
    /**
     * Load gene expression data downloaded from the data portal.
     * @param fileName
     * @return
     * @throws IOException
     */
    public Map<String, Map<String, Double>> loadGeneExpFromPortal(String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = fu.readLine();
        Map<String, Map<String, Double>> geneToExp = new HashMap<String, Map<String,Double>>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String sample = tokens[0].substring(0, 12);
            String gene = tokens[1];
            String value = tokens[2];
            double exp = MathUtilities.log2(new Double(value));
            Map<String, Double> map = geneToExp.get(gene);
            if (map == null) {
                map = new HashMap<String, Double>();
                geneToExp.put(gene, map);
            }
            map.put(sample, exp);
        }
        fu.close();
        return geneToExp;
    }
    
    public void zTransformBasedOnSample(Map<String, Map<String, Double>> geneToExp) {
        Map<String, SummaryStatistics> sampleToStat = new HashMap<String, SummaryStatistics>();
        for (String gene : geneToExp.keySet()) {
            Map<String, Double> sampleToExp = geneToExp.get(gene);
            for (String sample : sampleToExp.keySet()) {
                Double value = sampleToExp.get(sample);
                SummaryStatistics stat = sampleToStat.get(sample);
                if (stat == null) {
                    stat = new SummaryStatistics();
                    sampleToStat.put(sample, stat);
                }
                stat.addValue(value);
            }
        }
        for (String gene : geneToExp.keySet()) {
            Map<String, Double> sampleToExp = geneToExp.get(gene);
            for (String sample : sampleToExp.keySet()) {
                Double value = sampleToExp.get(sample);
                SummaryStatistics stat = sampleToStat.get(sample);
                Double transformed = (value - stat.getMean()) / stat.getStandardDeviation();
                sampleToExp.put(sample, transformed);
            }
        }
    }
    
}
