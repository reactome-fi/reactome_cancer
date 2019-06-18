/*
 * Created on Apr 23, 2012
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.funcInt.Protein;
import org.reactome.r3.fi.SurvivalAnalysisHelper;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to analyze data sets downloaded from TCGA
 * @author gwu
 *
 */
public class ICGCPancreaticCancerAnalyzer {
    private final String DIR_NAME = R3Constants.DATA_SET_DIR + "ICGC/";
    private FileUtility fu = new FileUtility();
    
    public ICGCPancreaticCancerAnalyzer() {
    }
    
    @Test
    public void doSurvivalAnalysis() throws Exception {
        Map<String, Set<String>> sampleToMutatedGenes = loadSampleToMutatedGenes();
        System.out.println("Total samples: " + sampleToMutatedGenes.size());
        Set<String> allGenes = InteractionUtilities.grepAllGenes(sampleToMutatedGenes);
        System.out.println("Total mutated genes: " + allGenes.size());
        String moduleGenesFile = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + "data_published_in_nature/GenesInModules6And7.txt";
        Set<String> moduleGenes = fu.loadInteractions(moduleGenesFile);
        System.out.println("Merged OV Exome module: " + moduleGenes.size());
        int mutatedSamples = 0;
        for (String sample : sampleToMutatedGenes.keySet()) {
            Set<String> mutatedGenes = sampleToMutatedGenes.get(sample);
            Set<String> shared = InteractionUtilities.getShared(mutatedGenes, moduleGenes);
            if (shared.size() > 0)
                mutatedSamples ++;
        }
        System.out.println("Total samples mutated in the OV exome module: " + mutatedSamples);
        
        TCGAOvarianCancerAnalyzer ovAnalyzer = new TCGAOvarianCancerAnalyzer();
        SurvivalAnalysisHelper survivalHelper = CancerAnalysisUtilitites.getSurvivalAnalysisHelper();
        String clinFileName = DIR_NAME + "Pancreatic_Cancer-QCMG-AU/SimpleClinInfo.txt";
        ovAnalyzer.doSurvivalAnalysisForMutationModule(survivalHelper,
                                                       sampleToMutatedGenes, 
                                                       moduleGenes, 
                                                       clinFileName, 
                                                       "Pancreatic_Cancer-QCMG-AU");
    }
    
    /**
     * Load sample to mutated genes.
     * @return
     * @throws IOException
     */
    private Map<String, Set<String>> loadSampleToMutatedGenes() throws IOException {
        String srcFileName = DIR_NAME + "Pancreatic_Cancer-QCMG-AU/snp_QCMGPancreasWithGenes.txt";
        String[] types = new String[] {
                "complex_indel",
                "essential_splice_site,intronic",
                "frameshift_coding",
                "frameshift_coding,splice_site",
                "splice_site,3prime_utr",
                "splice_site,5prime_utr",
                "stop_gained",
                "stop_gained,splice_site",
                "synonymous_coding"
        };
        Set<String> allowedTypes = new HashSet<String>(Arrays.asList(types));
        Map<String, Set<String>> sampleToGenes = new HashMap<String, Set<String>>();
        fu.setInput(srcFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 54)
                continue; // The gene column is 54
            if (!allowedTypes.contains(tokens[15].toLowerCase()))
                continue;
            InteractionUtilities.addElementToSet(sampleToGenes,
                                                 tokens[30],
                                                 tokens[53]);
        }
        return sampleToGenes;
    }
    
    @Test
    public void checkConsequenceTypes() throws IOException {
        String targetFileName = DIR_NAME + "Pancreatic_Cancer-QCMG-AU/snp_QCMGPancreasWithGenes.txt";
        Set<String> types = new HashSet<String>();
        fu.setInput(targetFileName);
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            types.add(tokens[15]);
        }
        fu.close();
        List<String> list = new ArrayList<String>(types);
        Collections.sort(list);
        for (String type : list)
            System.out.println(type);
    }
    
    /**
     * File snp_QCMGPancreatic.txt doesn't have a gene symbol column. Use this method to add
     * this column by mapping EMSEBML transcript id to uniprot to gene symbol.
     * @throws IOException
     */
    @Test
    public void addGeneSymbolColumnToSNPQCMGPancreas() throws IOException {
        // This file was generated by method UniProtAnalyzer.generateENSEBMLToUniMap()
        String fileName = R3Constants.RESULT_DIR + "ensembl2uni.txt";
        Map<String, Set<String>> ensemblToUniProt = fu.loadSetMap(fileName);
        Map<String, Protein> accessToProtein = new UniProtAnalyzer().generateUniAccToProteinMap();
        String srcFileName = DIR_NAME + "Pancreatic_Cancer-QCMG-AU/snp_QCMGPancreas.txt";
        String targetFileName = DIR_NAME + "Pancreatic_Cancer-QCMG-AU/snp_QCMGPancreasWithGenes.txt";
        fu.setInput(srcFileName);
        String line = fu.readLine();
        fu.setOutput(targetFileName);
        fu.printLine(line + "\tGene Symbol");
        StringBuilder builder = new StringBuilder();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String ensemblId = tokens[18];
            Set<String> uniProts = ensemblToUniProt.get(ensemblId);
            if (uniProts == null) {
                fu.printLine(line);
                continue;
            }
            builder.setLength(0);
            for (String id : uniProts) {
                Protein protein = accessToProtein.get(id);
                if (protein == null || protein.getShortName() == null)
                    continue;
                if (builder.length() > 0)
                    builder.append(",");
                builder.append(protein.getShortName());
            }
            if (builder.length() == 0)
                fu.printLine(line);
            else
                fu.printLine(line + "\t" + builder.toString());
        }
        fu.close();
    }
    
    /**
     * This method is used to create a simple survival file.
     * @throws IOException
     */
    @Test
    public void processSurvivalInfo() throws IOException {
//        String srcFileName = DIR_NAME + "Pancreatic_Cancer-QCMG-AU/sample_QCMGPancreas.txt";
//        String targetFileName = DIR_NAME + "Pancreatic_Cancer-QCMG-AU/SimpleClinInfo.txt";
        String srcFileName = DIR_NAME + "Pancreatic_Cancer-OICR-CA/sample_oicrPanc.txt";
        String targetFileName = DIR_NAME + "Pancreatic_Cancer-OICR-CA/SimpleClinInfo.txt";    
        StringBuilder builder = new StringBuilder();
        fu.setInput(srcFileName);
        fu.setOutput(targetFileName);
        String line = fu.readLine();
        fu.printLine("Sample\tOSEVENT\tOSDURATION");
        Set<String> samples = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (samples.contains(tokens[1]))
                continue;
            builder.append(tokens[1]).append("\t");
            if (tokens[2].equals("alive"))
                builder.append(0);
            else
                builder.append(1);
            builder.append("\t").append(tokens[8]);
            fu.printLine(builder.toString());
            builder.setLength(0);
            samples.add(tokens[1]);
        }
        fu.close();
    }
    
}
