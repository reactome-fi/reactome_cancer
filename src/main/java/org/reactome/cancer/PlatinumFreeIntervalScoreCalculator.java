/*
 * Created on Nov 9, 2009
 *
 */
package org.reactome.cancer;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.r3.graph.GraphComponent;
import org.reactome.r3.graph.GraphComponent.ScoreCalculator;
import org.reactome.r3.util.MathUtilities;

public class PlatinumFreeIntervalScoreCalculator implements ScoreCalculator {
    private Map<String, Double> sampleToPlatinum;
    private Map<String, Map<String, Double>> geneToExp;
    private boolean isForNegative = false;
    // To control is sample id should be processed
    private boolean useSampleIdDirectly;
    
    public PlatinumFreeIntervalScoreCalculator() {
    }
    
    public void setUseSampleIdDirectly(boolean flag) {
        this.useSampleIdDirectly = flag;
    }
    
    
    public Map<String, Double> getSampleToPlatinum() {
        return sampleToPlatinum;
    }



    public void setSampleToPlatinum(Map<String, Double> sampleToPlatinum) {
        this.sampleToPlatinum = sampleToPlatinum;
    }



    public Map<String, Map<String, Double>> getGeneToExp() {
        return geneToExp;
    }



    public void setGeneToExp(Map<String, Map<String, Double>> geneToExp) {
        this.geneToExp = geneToExp;
    }

    public void setIsForNegative(boolean flag) {
        this.isForNegative = flag;
    }

    public double calculateScore(GraphComponent comp) {
        Set<String> genes = comp.getAllNodes();
        return calculateScore(genes);
    }
    
    public double calculateScore(Collection<String> genes) {
     // Need to get average
        Map<String, Double> sampleToAvg = calculateAverage(genes);
        double score = calculateCorBetweenExpAndPlatInverval(sampleToAvg, 
                                                             sampleToPlatinum);
        if (isForNegative)
            score = -score;
        return score;
    }
    
    public Map<String, Double> calculateScoreForGene() {
        Map<String, Double> geneToScore = new HashMap<String, Double>();
        for (String gene : geneToExp.keySet()) {
            Map<String, Double> sampleToExp = geneToExp.get(gene);
            Double score = calculateCorBetweenExpAndPlatInverval(sampleToExp, sampleToPlatinum);
            geneToScore.put(gene, score);
        }
        return geneToScore;
    }
    
    protected double calculateCorBetweenExpAndPlatInverval(Map<String, Double> sampleToExp,
                                                           Map<String, Double> sampleToPlatinum) {
        List<Double> platInterval = new ArrayList<Double>();
        List<Double> expValues = new ArrayList<Double>();
        for (String sample : sampleToPlatinum.keySet()) {
            Double platInt = sampleToPlatinum.get(sample);
            if (hasExpressionValue(sample, sampleToExp)) {
                double expValue = getExpressionValue(sample, sampleToExp);
                platInterval.add(platInt);
                expValues.add(expValue);
            }
        }
        return MathUtilities.calculatePearsonCorrelation(platInterval, expValues);
    }
    
    protected double getExpressionValue(String sample, Map<String, Double> sampleToValue) {
        if (useSampleIdDirectly) {
            return sampleToValue.get(sample);
        }
        for (String expSample : sampleToValue.keySet()) {
            if (expSample.startsWith(sample)) {
                if (expSample.length() == sample.length())
                    return sampleToValue.get(expSample);
                // Have to make sure the sample type is a cancer sample (sampleId < 10)
                String sub = expSample.substring(sample.length());
                int sampleId = Integer.parseInt(sub.substring(1, 3)); // get code 01 in -01A- 
                if (sampleId < 10)
                    return sampleToValue.get(expSample);
            }
        }
        return Double.NaN;
    }
    
    protected boolean hasExpressionValue(String sample,
                                       Map<String, Double> sampleToExp) {
        if (useSampleIdDirectly) {
            return sampleToExp.containsKey(sample);
        }
        for (String expSample : sampleToExp.keySet()) {
            if (expSample.startsWith(sample)) {
                if (sample.length() == expSample.length())
                    return true; // Has been pre-processed
                // Have to make sure the sample type is a cancer sample (sampleId < 10)
                String sub = expSample.substring(sample.length());
                int sampleId = Integer.parseInt(sub.substring(1, 3)); // get code 01 in -01A- 
                if (sampleId < 10)
                    return true;
            }
        }
        return false;
    }
    
    protected Map<String, Double> calculateAverage(Collection<String> genes) {
        Map<String, Double> sampleToAvg = new HashMap<String, Double>();
        // In case there is null
        int size = 0;
        for (String gene : genes) {
            Map<String, Double> sampleToExp = geneToExp.get(gene);
            if (sampleToExp == null)
                continue;
            size ++;
            for (String sample : sampleToExp.keySet()) {
                Double value = sampleToExp.get(sample);
                Double tmp = sampleToAvg.get(sample);
                if (tmp == null)
                    sampleToAvg.put(sample, value);
                else
                    sampleToAvg.put(sample, tmp + value);
            }
        }
        // Calculate average
        for (String sample : sampleToAvg.keySet()) {
            Double tmp = sampleToAvg.get(sample);
            sampleToAvg.put(sample, tmp / size);
        }
        return sampleToAvg;
    }
    
}
