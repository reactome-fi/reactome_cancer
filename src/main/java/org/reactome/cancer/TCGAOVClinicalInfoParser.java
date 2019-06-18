/*
 * Created on Feb 2, 2010
 *
 */
package org.reactome.cancer;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.jdom.xpath.XPath;
import org.junit.Test;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to parse XML files for TCGA OV clinical information.
 * @author wgm
 *
 */
public class TCGAOVClinicalInfoParser {
    private DateFormat format = new SimpleDateFormat("MM/dd/yyyy");
    private Set<String> drugNames = new HashSet<String>();
    
    /**
     * Parse an input XML file to output file.
     * @param input
     * @param out
     * @throws Exception
     */
    public void parse(File file,
                      FileUtility fu) throws Exception {
        StringBuilder text = new StringBuilder();
        SAXBuilder builder = new SAXBuilder();
        Document doc = builder.build(file);
        Element root = doc.getRootElement();
        // Get the barcode for the patient
        String path = "PATIENT/BCRPATIENTBARCODE";
        Element barcode = (Element) XPath.selectSingleNode(root, path);
        text.append(barcode.getText()).append("\t");
        // Get the surgical time
        path = "PATIENT/SURGERIES/SURGERY";
        Element surgery = (Element) XPath.selectSingleNode(root, path);
        String year = surgery.getChildText("YEAROFPROCEDURE");
        String month = surgery.getChildText("MONTHOFPROCEDURE");
        String day = surgery.getChildText("DAYOFPROCEDURE");
        String date = createDateString(year, month, day);
        text.append(date).append("\t");
        // Get the last follow up
        String lastFollowup = getDateString(root, "LASTFOLLOWUP");
        text.append(lastFollowup).append("\t");
        // Get the date of progression
        String progression = getDateString(root, "TUMORPROGRESSION");
        text.append(progression == null ? "" : progression).append("\t");
        // Get the recurrence 
        String recurrence = getDateString(root, "TUMORRECURRENCE");
        text.append(recurrence == null ? "" : recurrence).append("\t");
        // Get the date for progree or recurrence
        Date relapseDate = null;
        if (progression != null)
            relapseDate = createDate(progression);
        else if (recurrence != null)
            relapseDate = createDate(recurrence);
        // Get the last time platinum 
        Date lastTimePlat = getLastTimePlat(root, relapseDate);
        text.append(lastTimePlat == null ? "" : format.format(lastTimePlat)).append("\t");
        // Check if there is relapse
        text.append(relapseDate == null ? "0" : "1").append("\t");
        if (lastTimePlat == null) {
            text.append("");
        }
        else {
            // Check plat-free-time
            double platFreeTime = getPlatFreeTime(lastTimePlat, 
                                                  relapseDate, 
                                                  createDate(lastFollowup));
            text.append(String.format("%.1f", platFreeTime));
        }
        fu.printLine(text.toString());
    }
    
    private double getPlatFreeTime(Date lastTimePlat, 
                                   Date replapseDate,
                                   Date lastFollowupDate) {
        long diff = 0;
        if (replapseDate != null)
            diff = replapseDate.getTime() - lastTimePlat.getTime();
        else
            diff = lastFollowupDate.getTime() - lastTimePlat.getTime();
        // Diff in ms
        double rtn = diff / (1000.0 * 60 * 60 * 24);
        return rtn  * 12.0 / 365.0;
    }
    
    private Date getLastTimePlat(Element root,
                                 Date relapseDate) throws Exception {
        Element drugsElm = (Element) XPath.selectSingleNode(root, "PATIENT/DRUGS");
        // Handle platinum drug only
        List list = drugsElm.getChildren("DRUG");
        Date rtnDate = null;
        for (Iterator it = list.iterator(); it.hasNext();) {
            Element drugElm = (Element) it.next();
            String drugName = drugElm.getChildText("DRUGNAME");
            drugNames.add(drugName);
            if (!drugName.endsWith("platin"))
                continue;
            String day = drugElm.getChildText("DAYOFDRUGTREATMENTEND");
            String month = drugElm.getChildText("MONTHOFDRUGTREATMENTEND");
            String year = drugElm.getChildText("YEAROFDRUGTREATMENTEND");
            if (day.length() == 0 || month.length() == 0 || year.length() == 0)
                continue;
            String dateString = createDateString(year, month, day);
            // Check if this date should be used
            Date platDate = createDate(dateString);
            if (relapseDate == null)
                rtnDate = platDate;
            else if (platDate.before(relapseDate)) {
                rtnDate = platDate;
            }
        }
        return rtnDate;
    }
    
    private void printHeaders(FileUtility fu) throws IOException {
        fu.printLine("BCRPATIENTBARCODE\t1st_SURGERIES/SURGERY/DATEOFPROCEDURE" +
        		"\tDATEOFLASTFOLLOWUP\tDATEOFTUMORPROGRESSION\tDATEOFTUMORRECURRENCE" +
        		"\tLAST_DATE_PRIMARY_PLATINUM_TREATMENT\tPROGRESSION_STATUS\tPLATINUM_FREE_INTERVAL_MONTHS");
    }
    
    private String getDateString(Element root,
                                 String name) throws JDOMException {
        String path = "PATIENT/DAYOF" + name;
        String day = ((Element)XPath.selectSingleNode(root, path)).getText();
        path = "PATIENT/MONTHOF" + name;
        String month = ((Element)XPath.selectSingleNode(root, path)).getText();
        path = "PATIENT/YEAROF" + name;
        String year = ((Element)XPath.selectSingleNode(root, path)).getText();
        if (day.length() == 0 || month.length() == 0 || year.length() == 0)
            return null;
        return createDateString(year, month, day);
    }
    
    private Date createDate(String date) throws ParseException {
        return format.parse(date);
    }
    
    private String createDateString(String year, String month, String day) {
        return month + "/" + day + "/" + year;
    }
    
    
    @Test
    public void testParse() throws Exception {
        String dirName = TCGAOvarianCancerAnalyzer.OVARIAN_DIR_NAME + "ClinicalInfo012110/";
        List<File> files = new ArrayList<File>();
        String[] subDirNames = new String[] {
                "intgen.org_OV.bio.Level_1.19.5.0",
                "intgen.org_OV.bio.Level_1.18.6.0",
                "intgen.org_OV.bio.Level_1.17.7.0"
        };
        for (String name : subDirNames) {
            File dir = new File(dirName + name);
            for (File file : dir.listFiles()) {
                if (file.getName().endsWith(".xml"))
                    files.add(file);
            }
        }
        FileUtility fu = new FileUtility();
        String output = dirName + "ClinBatch17_19.txt";
        fu.setOutput(output);
        printHeaders(fu);
        for (File input : files) {
            System.out.println("File: " + input.getName());
            parse(input, fu);
        }
        fu.close();
        // Check drug names
        for (String drug : drugNames)
            System.out.println(drug);
    }
}
