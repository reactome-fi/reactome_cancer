/*
 * Created on Feb 10, 2012
 *
 */
package org.reactome.cancer;

import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

import org.junit.Test;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.GeneExpressionDataSet;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to do some data analysis for Pat Shaw's group around OV gene expression data sets.
 * @author gwu
 *
 */
public class PatShawDataAnalyzer {
    private final String dirName = R3Constants.DATA_SET_DIR + "pat_shaw/";
    private FileUtility fu = new FileUtility();
    
    public PatShawDataAnalyzer() {
        
    }
    
    /**
     * Process clinical informaiton file.
     * @throws Exception
     */
    @Test
    public void processClinicaInfo() throws Exception {
        String srcFile = dirName + "clinical data profiled cases without MRN-2.txt";
        String outFileName = dirName + "Clinica_OV.txt";
        fu.setInput(srcFile);
        fu.setOutput(outFileName);
        String line = fu.readLine();
        fu.printLine("Sample\tOSEVENT\tOSDURATION");
        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd");
        long millInDay = (24 * 60 * 60 * 1000);
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            // Date of diagnosis
            Date diagnosisDate = dateFormat.parse(tokens[1]);
            // Date of Most Recent Contact
            Date mostRecentDate = dateFormat.parse(tokens[15]);
            // Calculate OS
            long timeDiff = mostRecentDate.getTime() - diagnosisDate.getTime();
            int days = (int) (timeDiff / millInDay);
            double months = (days / 30.0d);
            // Patient Status
            String status = tokens[13];
            String statusValue = "";
            if (status.startsWith("Dead"))
                statusValue = "1";
            else if (status.startsWith("Alive"))
                statusValue = "0";
            fu.printLine(tokens[0] + "\t" +
                         statusValue + "\t" +
                         String.format("%.2f", months));
        }
        fu.close();
    }
    
    @Test
    public void generateSampleToGeneExpClusters() throws IOException {
        String clusterFileName = CancerResequenceDataSetAnalyzer.OVARIAN_DIR_NAME + "MCL_Clusters_FIsWithGeneExpAbsCorr_102510_I50.txt";
        String fiToCorFileName = CancerResequenceDataSetAnalyzer.OVARIAN_DIR_NAME + "FIsWithGeneExpAbsCorr_100510.txt";
        
//        String geneExpFileName = dirName + "Milea_Shaw_43Cases_HGSC_log_z_transformed_021412.txt";
//        String output = dirName + "Milea_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_log_z_transformed_021412.txt";
        
//        String geneExpFileName = dirName + "Milea_Shaw_43Cases_HGSC_z_transformed_021012.txt";
//        String output = dirName + "Milea_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_z_transformed_021012.txt";
        
        String geneExpFileName = dirName + "Milea_Shaw_samples_021412.txt";
        String output = dirName + "Milea_SampleToGeneExpMCL_Clusters_FIsWithGeneExpAbsCorr_021412.txt";
        
        CancerGeneExpressionCommon helper = new CancerGeneExpressionCommon();
        Map<String, Map<String, Double>> geneToSampleToValue = helper.loadGeneExp(geneExpFileName);
        helper.generateSampleToGeneExpClusters(clusterFileName,
                                               fiToCorFileName,
                                               geneToSampleToValue,
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_SIZE_CUTOFF,
                                               OvarianMCLGeneExpModuleAnalyzer.MCL_MEAN_CUTOFF,
                                               output);
    }
    
    @Test
    public void processMileaShaw43CasesHGSCFile() throws Exception {
        Map<String, String> idToGene = loadProbesetToGene();
        System.out.println("Total mappable ids: " + idToGene.size());
//        String srcFileName = dirName + "Milea_Shaw_43Cases_HGSC.txt";
//        String destFileName = dirName + "Milea_Shaw_43Cases_HGSC_log_z_transformed_021412.txt";
        
        String srcFileName = dirName + "Milea_Shaw_samples.txt";
        String destFileName = dirName + "Milea_Shaw_samples_021412.txt";
        
        GSEMatrixDataHandler matrixDataHandler = new GSEMatrixDataHandler();
        GeneExpressionDataSet expressionDataset = matrixDataHandler.loadGeneExpressionDataSet(srcFileName,
                                                                                              54615, // All probesets starting with AFFX have been excluded
                                                                                              1);
        System.out.println("Total features: " + expressionDataset.getFeatureList().size());
        System.out.println("Total samples: " + expressionDataset.getSampleList().size());
        matrixDataHandler.mapProbesetToGenes(expressionDataset, idToGene);
        // Do a log-transformation first
//        expressionDataset.logTransformation();
        expressionDataset = matrixDataHandler.averageValuesForSameGenes(expressionDataset);
        // zscore transformation
//        expressionDataset.logTransformation();
//        expressionDataset.zscoreTansformation();
        expressionDataset.export(destFileName);
    }
    
    private Map<String, String> loadProbesetToGene() throws IOException {
        String srcFileName = dirName + "Normalized-all-George_Shaw_data_all.txt";
        fu.setInput(srcFileName);
        String line = fu.readLine();
        Map<String, String> idToGene = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 2)
                continue; // Some funny lines
            if (tokens[1].length() == 0)
                continue; // Cannot be mapped
            if (tokens[1].contains("///"))
                continue; // Make sure one probeId to one gene only
            idToGene.put(tokens[0], tokens[1]);
        }
        fu.close();
        return idToGene;
    }
    
}
