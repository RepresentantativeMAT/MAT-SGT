/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package measure;


import br.ufsc.model.AttributeValue;
import br.ufsc.model.Centroid;
import br.ufsc.model.MultipleAspectTrajectory;
import br.ufsc.model.Point;
import br.ufsc.model.STI;
import br.ufsc.model.SemanticAspect;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author vanes
 */
public class SimilarityMeasure {

    private Map<Object, Float> weights;
    private Map<SemanticAspect, Float> thresholds;
//    private String patternTemporalData;

    private double recallRT_T = 0;

    public SimilarityMeasure() {
        this.weights = new HashMap<Object, Float>();
        this.thresholds = new HashMap<SemanticAspect, Float>();
    }


    public void clear() {
        this.weights.clear();
        this.thresholds.clear();
    }

    public void setWeight(Object attribute, float weight) {
            this.weights.put(attribute, weight);
    }

    public double getWeight(Object attribute) {
        try {
            if(attribute instanceof STI){
//                System.err.println("Deu certo o weight para Time: "+this.weights.get("TIME"));
                return this.weights.get("TIME");
            }else{
//                System.err.println("weight of: "+attribute+" = "+this.weights.get(attribute));
                return this.weights.get(attribute);
            }
        } catch (Exception e) {
            System.err.println("Error in getWeight for feature: '" + attribute + "' (weights: " + this.weights + ")");
            throw new NullPointerException();
        }
    }
    
    
    public float getAllWeight() {
        float sumWeight = 0.0f;
        

        for (Map.Entry<Object, Float> eachWeight : weights.entrySet()) {
            sumWeight += eachWeight.getValue();
            
        }
        return sumWeight;
    }

    public void setThreshold(SemanticAspect att, float threshold) {
        this.thresholds.put(att, threshold);
//        System.out.println("Lista de thresholds");
    }

    public double getThreshold(SemanticAspect att) {
        if (thresholds.isEmpty()) {
            System.err.println("threshold list is Empty");
        }

        try {
            return this.thresholds.get(att);
        } catch (Exception e) {
            System.err.println("Error in getThreshold for attribute: '" + att.getName() + "' (thresholds: " + this.thresholds + ")");
            throw new NullPointerException();
        }
    }

    public double getRecallRT_T() {
        return recallRT_T;
    }

    
    /**
     * Compute how much RT is capturing on the information contained in each original Trajectory (T). This value is given by parity between T and RT divided by the size of the T.
     * 
     * @param RT - representative trajectory 
     * @param t2 - an input trajectory
     * @return the recall value
     * @throws ParseException 
     */
    public double recallOf(MultipleAspectTrajectory RT, MultipleAspectTrajectory t2) throws ParseException {
        recallRT_T = 0;
        //Matrix RT x T
        double[][] scores = new double[RT.getPointList().size()][t2.getPointList().size()];

        //for each point of the input trajectory, analyse with all points of the RT
        for (int j = 0; j < t2.getPointList().size(); j++) {
            double maxCol = 0;

            for (int i = 0; i < RT.getPointList().size(); i++) {
                scores[i][j] = this.score((Centroid) RT.getPointList().get(i), t2.getPointList().get(j));
                maxCol = scores[i][j] > maxCol ? scores[i][j] : maxCol; //identify the point of the RT more similar (max score) with the actual point of input T.
            }

            recallRT_T += maxCol;
        }
        return (recallRT_T) / (t2.getPointList().size());

    }

    /**
     * Compute score value between two points - one of RT and one of T2
     * the score value is the sum of each aspect match in context. 
     * the aspect match value is considered following:
     *  - spatial aspect - when the point of T2 is into the same cell of RT centroid.
     *  - temporal aspect - when the temporal value of T2 is into some interval in RT or if it is equal an ponctual occurence.
     *  - semantic dimension - categorical aspects: when the value of T2 is into some interval of RT.
     *  - semantic dimension - numerical aspects: when the value of T2 is into the defined threshold of difference values with RT value.
     * 
     * @param p1 - Centroid point
     * @param p2 - a point of T2
     * @return score value
     * @throws ParseException 
     */
    private double score(Centroid p1, Point p2) throws ParseException {
        double score = 0;
        //if centroid contains p of T, then the match occurs
        if (p1.getPointListSource().contains(p2)) {
            score += (getWeight("SPATIAL"));
        }

        matchTemporal:
        {
            double match = 0;
            //if the time occurrence of p of T is into the centroid interval, the match occurs
            if (p1.getSti().getInterval().isInInterval(p2.getTime().getStartTime())){
                match = p1.getSti().getProportion();
            }

            score += match * getWeight(p1.getSti());
            
            
        }

        for (AttributeValue atvP1 : p1.getListAttrValues()) { //loop on list of attribute values of centroid
            AttributeValue tempAttP2 = p2.getAttributeValue(atvP1.getAttibute()) != null ? p2.getAttributeValue(atvP1.getAttibute()) : null;
            double tempSemanticMatch = computeMatch(atvP1, tempAttP2);
            score += tempSemanticMatch;

        }
        return score;
    }

    /**
     * Compute the match value between two attribute values (RT and T2)
     * the match value is the final value of each attribute - match value * defined weight
     * 
     * @param rep
     * @param atv
     * @return the value of the match
     */
    public double computeMatch(AttributeValue rep, AttributeValue atv) {
        double match = 0;
        if(atv == null || rep == null)
            return 0;


                if (rep.getValue() instanceof Map) {
                // case of semantic - categorical
                HashMap<String, Double> valuesRT = (HashMap) rep.getValue();
                if (valuesRT.containsKey(((String) atv.getValue()).toUpperCase())) {
                    match = valuesRT.get(((String) atv.getValue()).toUpperCase());
                }
            } else { // numerical values
                
                if(rep.getNumericalValueSD() != 0)
                    match = Math.abs(Double.parseDouble((String)rep.getValue()) - Double.parseDouble((String)atv.getValue())) <= (rep.getNumericalValueSD() * 2.5) ? 1.0 : 0;
                else //default value
                    match = Math.abs(Double.parseDouble((String)rep.getValue()) - Double.parseDouble((String)atv.getValue())) <= 10 ? 1.0 : 0;
//                
//                System.out.println("Numerical value: "+atv.getValue()+" -- RT value: "+rep.getValue()+" --Threshold value:  "+rep.getNumericalValueSD()*2.5+" -- Match value: "+match);
            }

        return match * getWeight(rep.getAttibute());
    }


}
