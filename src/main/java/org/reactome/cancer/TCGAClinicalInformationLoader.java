/*
 * Created on May 11, 2010
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.reactome.r3.util.FileUtility;

/**
 * This class is used to load TCGA clinical information.
 * @author wgm
 *
 */
public class TCGAClinicalInformationLoader {
    private FileUtility fu;
    
    public TCGAClinicalInformationLoader() {
        fu = new FileUtility();
    }
    
    public Map<String, TCGAClinicalInfo> loadClinicalInfo(String fileName) throws IOException {
        Map<String, TCGAClinicalInfo> rtn = new HashMap<String, TCGAClinicalInfo>();
        Map<String, Double> sampleToRate = new HashMap<String, Double>();
        fu.setInput(fileName);
        String line = fu.readLine();
        // Get the index
        int eventIndex = 0;
        int rateIndex = 0;
        String[] tokens = line.split("\t");
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals("OSEVENT"))
                eventIndex = i;
            else if (tokens[i].equals("OSDURATION_MONTHS"))
                rateIndex = i;
        }
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            String sample = tokens[0];
            TCGAClinicalInfo info = new TCGAClinicalInfo();
            info.setSample(sample);
            rtn.put(sample, info);
            if (tokens.length < rateIndex + 1)
                continue;
            String interval = tokens[rateIndex];
            if (interval.length() > 0)
                info.setOsDuration(new Double(interval));
            String type = tokens[eventIndex];
            if (type.length() > 0) {
                if (type.equals("0"))
                    info.setOsEvent(Boolean.FALSE);
                else if (type.equals("1"))
                    info.setOsEvent(Boolean.TRUE);
            }
        }
        fu.close();
        return rtn;
    }
    
    public Map<String, Double> loadSampleToSurvivalRate(boolean excludeLastFlowup,
                                                        String fileName) throws IOException {
        Map<String, Double> sampleToRate = new HashMap<String, Double>();
        fu.setInput(fileName);
        String line = fu.readLine();
        // Get the index
        int eventIndex = 0;
        int rateIndex = 0;
        String[] tokens = line.split("\t");
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals("OSEVENT"))
                eventIndex = i;
            else if (tokens[i].equals("OSDURATION_MONTHS"))
                rateIndex = i;
        }
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            String sample = tokens[0];
            if (tokens.length < rateIndex + 1)
                continue;
            if (!excludeLastFlowup) {
                String interval = tokens[rateIndex];
                if (interval.length() > 0) {
                    double value = Double.parseDouble(interval);
                    if (value >= 0) // Want to keep only meaningful values
                        sampleToRate.put(sample, value);
                }
            }
            else {
                String type = tokens[eventIndex];
                if (type.equals("1")) {
                    String interval = tokens[rateIndex];
                    if (interval.length() > 0) {
                        double value = Double.parseDouble(interval);
                        if (value > 0) // Want to keep only meaningful values
                            sampleToRate.put(sample, value);
                    }
                }
                else {
                    String interval = tokens[rateIndex];
                    if (interval.length() > 0) {
                        double value = Double.parseDouble(interval);
                        if (value > 40.0d) // Want to keep only meaningful values
                            sampleToRate.put(sample, value);
                    }
                }
            }
        }
        fu.close();
        return sampleToRate;
    }
    
    public Map<String, Double> loadSampleToPlatinumFreeInterval(boolean excludeLastFlowUp,
                                                                String fileName) throws IOException {
        Map<String, Double> sampleToInterval = new HashMap<String, Double>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String sample = tokens[0];
            // Second to last is PLATINUM_FREE_INTERVAL_MONTHS.
            // However index 46 should be used since some rows don't have enough values
            if (tokens.length < 46)
                continue;
            if (!excludeLastFlowUp) {
                String interval = tokens[45];
                if (interval.length() > 0) {
                    double value = Double.parseDouble(interval);
                    if (value >= 0) // Want to keep only meaningful values
                        sampleToInterval.put(sample, value);
                }
            }
            else {
                String type = tokens[44];
                if (type.equals("1")) {
                    String interval = tokens[45];
                    if (interval.length() > 0) {
                        double value = Double.parseDouble(interval);
                        if (value > 0) // Want to keep only meaningful values
                            sampleToInterval.put(sample, value);
                    }
                }
            }
        }
        fu.close();
        return sampleToInterval;
    }
    
    public Map<String, String> loadSampleToPrimaryTherapyOutcome(String fileName) throws IOException {
        Map<String, String> sampleToOutput = new HashMap<String, String>();
        fu.setInput(fileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String sample = tokens[0];
            if (tokens.length < 26)
                continue;
            String output = tokens[25];
            if (output.length() > 0) {
                sampleToOutput.put(sample, output);
            }
        }
        fu.close();
        return sampleToOutput;
    }
}
