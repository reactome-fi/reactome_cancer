/*
 * Created on Sep 13, 2006
 *
 */
package org.reactome.r3;

import java.io.IOException;

import junit.framework.TestCase;

import org.reactome.r3.util.FileUtility;

/**
 * Analyzer for R based results
 * @author guanming
 *
 */
public class RDataAnalyzer extends TestCase {
    private final String RESULT_DIR = "results/";
    
    public void logitTestResults() throws IOException {
        String fileName = RESULT_DIR + "PantherDataWithoutReactome091106_1.arff";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        double x, e, p;
        int relT = 0;
        int relF = 0;
        int trueT = 0;
        int falseT = 0;
        int trueF = 0;
        int falseF = 0;
        boolean isT = false;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("@"))
                continue;
            if (line.length() == 0)
                continue;
            String[] tokens = line.split(",");
            x = -2.938456d;
            if (tokens[0].equals("true")) {
                relT ++;
                isT = true;
            }
            else {
                relF ++;
                isT = false;
            }
            if (tokens[1].equals("true"))
                x += 1.828103d;
            if (tokens[2].equals("true"))
                x += 0.879601;
            if (tokens[3].equals("true"))
                x += 1.609404;
            if (tokens[4].equals("none"))
                x -= 0.076597;
            else if (tokens[4].equals("pos"))
                x += 0.193847;
            if (!tokens[5].equals("?"))
                x += Double.parseDouble(tokens[5]) * 0.624362;
            if (!tokens[6].equals("?"))
                x += Double.parseDouble(tokens[6]) * 0.377490;
            e = Math.exp(x);
            p = e / (e + 1);
            if (p > 0.5) {
                if (isT)
                    trueT ++;
                else
                    falseT ++;
            }
            else {
                if (isT)
                    falseF ++;
                else
                    trueF ++;
            }
        }
        fu.close();
        System.out.printf("positive %d; predicate true %d, false %d%n", relT, trueT, falseF);
        System.out.printf("negative %d; predicate true %d, false %d%n", relF, falseT, trueF);
    }
    
}
