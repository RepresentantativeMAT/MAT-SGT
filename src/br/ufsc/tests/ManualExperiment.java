/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package br.ufsc.tests;





import br.ufsc.methods.MATSG_TimeSeq_version3;
import br.ufsc.methods.MATSGT;
import java.io.IOException;
import java.text.ParseException;

/**
 *
 * @author vanes
 */
public class ManualExperiment {

    public static String filename;
    public static String extension;
    public static String dir;

    public static void main(String[] args) throws IOException, ParseException, CloneNotSupportedException {

        dir = "datasets/";
//        dir += "Gowala/";
        dir += "RE/";
        filename = "Running_Example_v5";
//        filename = args[0];
        extension = ".csv";

        
        //informando lista de att a ser forçados como categoricos, mesmo contendo números
        String[] lstCategoricalsPreDefined = {"price"};
        for (int i = 0; i < lstCategoricalsPreDefined.length; i++) {
            lstCategoricalsPreDefined[i] = lstCategoricalsPreDefined[i].toUpperCase();
        }
        
        String SEPARATOR = ",";
        
        String[] valuesNulls = {"-1"};

        String[] lstIgnoreColumns = null;
//        String[] lstIgnoreColumns = {"label"};
//        for (int i = 0; i < lstIgnoreColumns.length; i++) {
//            lstIgnoreColumns[i] = lstIgnoreColumns[i].toUpperCase();
//        }

        float threshold_rc = 0.25f;
        float threshold_rv = 0.25f;
//        String patternDate = "yyyy-MM-dd HH:mm:SS.SSS";
        String patternDateIn = "?"; //For minutes time (integer value) inform '?' character
        
        MATSGT method = new MATSGT();
//            method.notConsiderNulls();
            method.setFilenameFullDataset(filename+"_complete");
            method.execute(dir, filename, extension, lstCategoricalsPreDefined, SEPARATOR, valuesNulls, lstIgnoreColumns, patternDateIn, threshold_rc, threshold_rv);
        
        
        
    }

}
