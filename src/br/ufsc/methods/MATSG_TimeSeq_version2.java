/*
    version 2 test: -- update on 19/04/23
    
    verificar como está minha medida de similaridade, quais os pontos está analisando
    -- identificar problema maior...
     --- se espaço? usar distância euclideana ao invés da célula? 
            ----> R: ok -- update to 2x tau spatial


     --- se temporal? usar feature para associar ao espacial?



    --- semantico? Está de fato calculando somente os scores conforme parametrizado do Rank? tau_rv

     --- Se sim? 
            --- Está reajustando os valores restante em um limiar de 0-1? 
                -----> Ok --> reajustado em normalizeRankingValues
     
        
 */
package br.ufsc.methods;

import br.ufsc.model.AttributeValue;
import br.ufsc.model.Centroid;
import br.ufsc.model.MultipleAspectTrajectory;
import br.ufsc.model.Point;
import br.ufsc.model.STI;
import br.ufsc.model.SemanticAspect;
import br.ufsc.model.SemanticType;
import br.ufsc.model.TemporalAspect;
import br.ufsc.util.Util;
import br.ufsc.util.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import measure.MUITAS;
import measure.SimilarityMeasure;

/**
 *
 * @author Vanessa
 */
public class MATSG_TimeSeq_version2 {
    // Attributes 

    // ------------- to Spatial division -- Dataset file information
    private String filename; //Filename of the dataset
    private String directory;//Directory of the dataset
    private String extension; //Extension of the filename

    private String SEPARATOR;

    // --- Define initial index value to semantic attributes
    private int INDEX_SEMANTIC = 3;
    private boolean considerNulls = true;

    private final String nullValue = "Unknown".toUpperCase();
    private final int THRESHOLD_TIME = 100;

    // --------------- to determine categoricals pre-defined values
    List<String> lstCategoricalsPD;
    List<String> lstIgnoreCols = null;
    String[] valuesNulls; //witch values are considered null in input dataset?

    //-- parameters to MAT-SG (defined by the user)
    private float threshold_rc; //To define relevant cells 
    private float threshold_rv; //To define relevant values in rank values, which values in rank are representative
    private float rc; //To define relevant values in rank values, which values in rank are representative

    // -- Load
    // For loading information from the dataset
    private static List<SemanticAspect> attributes; //List of all diferent attributes found in the dataset
    private static List<Point> points; //Points to be analysed
    private static Map<String, BitSet> spatialCellGrid; //Spatial grid array
    private static List<Point> pointsInCell; //List of all diferent points found in the each cell

    // format of input data
    private static SimpleDateFormat formatDate;
    // format of output data
    private static SimpleDateFormat formatDateOut;

// --------------------- AUX ----------------------
    private static int rId;
    private static String auxTid;

    // To model trajectory data
    private static MultipleAspectTrajectory trajectory; //Contain all points of a MAT
    private static List<MultipleAspectTrajectory> listTrajectories; //List of all MATs in the dataset
    private static MultipleAspectTrajectory betterRT; //Better Summarized MAT

    /// ----- Spatial Segmentation
    // To create the Spatial division    
    private static float spatialThreshold; //Maximum possible size for a cell
    private static double cellSizeSpace; //Size of each cell
    private float auxMaxZ; //

    /// ---- Summarization step
    // To create the Temporal summarization
    private ArrayList<Date> listTimesInCell; //List of all time marks in a cell -- update: option times in Date

    //aux representative MAT for ordenate
    private static List<Centroid> listRepPoint;

    private static MultipleAspectTrajectory representativeTrajectory; //Summarized MAT

    //aux to know the cell of each rp
    private String presentCell;

    // To create the Spatial summarization
    private double avgX, avgY;

    // To provide and compute runtime
    private Date initialTemp;

    private boolean dailyInfo;

    //Sum of each type of numerical attribute
    private static Map<String, List<Double>> sematicNumericFusionVal;
    //Sum of ocorrunces os each categorical attribute
    private static Map<Object, Map<String, Integer>> sematicCategoricalSummarizationVal;

    // For validation
    private String filenameFullDataset; //Filename of the dataset
    private static List<MultipleAspectTrajectory> listAllTrajectories; //List of all MATs in the dataset

    //pattern to number
    DecimalFormat formatNumber = new DecimalFormat("0.##", DecimalFormatSymbols.getInstance(Locale.US));

    /**
     * Reads the dataset file and creates the all the MATs
     *
     * @throws IOException
     */
    private void load() throws IOException, ParseException {

        java.io.Reader input = new FileReader(directory + filename + extension);
        BufferedReader reader = new BufferedReader(input);

        String datasetRow = reader.readLine();
        //To Get the header of dataset
        String[] datasetColumns = datasetRow.split(SEPARATOR);

        //To add all types of semantic attributes in the dataset, specified in the first line
        int order = 0;
        for (String s : Arrays.copyOfRange(datasetColumns, INDEX_SEMANTIC, datasetColumns.length)) {
            //when attr do not need to be ignored
            if (lstIgnoreCols == null || !lstIgnoreCols.contains(s.toUpperCase().trim())) {
                if (lstCategoricalsPD.contains(s.toUpperCase())) //If attribute was predefined as categorical
                {
                    attributes.add(new SemanticAspect(s.toUpperCase().trim(), order++, SemanticType.CATEGORICAL));
                } else {
                    attributes.add(new SemanticAspect(s.toUpperCase().trim(), order++));

                }
            } else {
                order++; //to skip column when it need to be ignored
            }

        }

        datasetRow = reader.readLine();

        //EoF - To get the trajectory data of dataset of each line
        while (datasetRow != null) {
//            System.out.println("linha: "+datasetRow);
            datasetColumns = datasetRow.toUpperCase().split(SEPARATOR);
            addAttributeValues(datasetColumns);
            datasetRow = reader.readLine();
        }
//        System.out.println("MATS: ");
//        for (MultipleAspectTrajectory mat: listTrajectories) {
//            System.out.println(mat);
//        }
//        System.exit(0);

        reader.close();

    }

    /**
     *
     * @param attrValues
     * @throws ParseException
     */
    private void addAttributeValues(String[] attrValues) throws ParseException {

        ++rId; //Id given to each data point 

        //Defines the semantic dimension as all attributes in predefined index to the end of line
        String[] semantics = Arrays.copyOfRange(attrValues, INDEX_SEMANTIC, attrValues.length);

        //All trajectory point follow the pattern:
        //id trajectory, coordinates (lat long), time, all semantic dimensions...
        // Follow the pattern add each MAT point in relative MAT
        if (formatDate != null) {
            addTrajectoryData(attrValues[0], attrValues[1].split(" "), formatDate.parse(attrValues[2]), semantics);
        } else {
            addTrajectoryData(attrValues[0], attrValues[1].split(" "), Util.convertMinutesToDate(Integer.parseInt(attrValues[2])), semantics);
        }

    } // end of addAttributeValue method

    /**
     * Add each MAT point in relative MAT object -- mapping input data to the
     * model predefined following O.O.
     *
     * @param tId - Id of MAT
     * @param coordinates - coordinates of point
     * @param time - time date of point
     * @param semantics - semantics attributes of point
     */
    private void addTrajectoryData(String tId, String[] coordinates, Date time, String[] semantics) {

        if (!tId.equals(auxTid)) { //IF the MAT is not created
            auxTid = tId;
            listTrajectories.add(new MultipleAspectTrajectory(Integer.parseInt(tId))); //Adds (Create) the new trajectory
            trajectory = listTrajectories.get(listTrajectories.size() - 1);
        }

        // aux values
        ArrayList<AttributeValue> attrs = new ArrayList<>();
        int ord = 0;
        SemanticAspect a;

        //Organizes the point semantic attributes
        for (String val : semantics) {

            a = findAttributeForOrder(ord++);
            if (a != null) { // This one will be NULL whether the columns is setted as ignored
//                System.out.println("Attribute order: " + a + " --> val: " + val);
//                if (Arrays.asList(valuesNulls).contains(val)) { // Define as Unknown null values
////                    System.out.println("Value null");
//                    val = nullValue;
//                } else 
                if (a.getType() != null
                        && a.getType().equals(SemanticType.CATEGORICAL)) { //if it is predefined as Categorical

                    val = "*" + Integer.parseInt(val); // Use character '*' to force the number value to be a categorical value
                }

                attrs.add(new AttributeValue(val.toUpperCase(), a));

            }

        }
        a = null; //clean memory 

        //Vanessa: parei aqui
        //Adds the MAT point to current MAT
        trajectory.addPoint(new Point(rId,// Point ID -- a sequencial value starting in 1...
                Double.parseDouble(coordinates[0]),
                Double.parseDouble(coordinates[1]),
                time,
                attrs));

        //Adds current MAT point to list of points
        points.add(trajectory.getLastPoint());

    } //end addTrajectoryData method

    ///------
    //1st step of MATSG --> Spatial Segmentation
    /**
     * add each point in the relative grid cell
     *
     * @param coordinates Coordinates of each trajectory point
     *
     */
    private static void allocateInSpaceCell(Point p) {

        //Get x,y of the point in the spatial grid
        String key = getCellPosition(p.getX(), p.getY());

        //Get id of the spatial grid cell
        BitSet rIds = spatialCellGrid.get(key);

        //If the cell doesn't exist
        if (rIds == null) {
            //Creates the cell and adds to the spatial grid
            rIds = new BitSet();
            rIds.set(p.getrId());
            spatialCellGrid.put(key, rIds);
        } else {
            //update on cell (spatial grid) adding this point
            rIds.set(p.getrId());
            spatialCellGrid.replace(key, rIds);
        }
    } //end allocateInSpaceCell method

    /**
     * Compute the cell position based on x and y divided by the cell size
     * predefined (cellSizeSpace)
     *
     * It calculates the cell position by: 1. dividing the x-coordinate by
     * cellSizeSpace and taking the floor value. The result is cast to an
     * integer: (int) Math.floor(x / cellSizeSpace). This determines the row
     * index of the cell position.
     *
     * 2. on the cell position, the column index is determined dividing the
     * y-coordinate by cellSizeSpace, taking the floor value, and casting it to
     * an integer: (int) Math.floor(y / cellSizeSpace).
     *
     * The floor is used to round down the division result to the nearest
     * integer value.
     *
     * @param x (x-coordinate)
     * @param y (y-coordinate)
     * @return Cell Position
     */
    private static String getCellPosition(double x, double y) {

        return ((int) Math.floor(x / cellSizeSpace)) + "," + ((int) Math.floor(y / cellSizeSpace));

    } // end getCellPosition method

    ///--------
    /// 2nd MATGS TimeSeq step --> Summarization
    //OK
    /**
     * Analyse all valid cell and for each cell -- identify the points in cell
     * (add all object Point) -- identify and list all times, i.e. the temporal
     * value of each point -- then summarize temporal aspect (by
     * summarizeTemporalAspect method)
     */
    public void identifyTimesInCell() {

        //Create iterator object of all spatial grid cells
        for (Map.Entry<String, BitSet> entry : spatialCellGrid.entrySet()) {
            String cellAnalyzed = entry.getKey();

            //Gets amount of points in the current cell
            int qntPoints = spatialCellGrid.get(cellAnalyzed).cardinality();
//            System.out.println("Cell Analyzed: " + cellAnalyzed + " | points: " + qntPoints + " | threshold RC: " + threshold_rc);
            if (qntPoints >= threshold_rc) { // IF number is at least a threshold RC
                resetValuesToSummarization();

                // Loop in all points of the cell
                for (int pointId = spatialCellGrid.get(cellAnalyzed).nextSetBit(0);
                        pointId >= 0;
                        pointId = spatialCellGrid.get(cellAnalyzed).nextSetBit(pointId + 1)) {
                    pointsInCell.add(points.get(pointId - 1));

                    //Temporal data 
                    listTimesInCell.add(points.get(pointId - 1).getTime().getStartTime()); // update: add start time (in Date) of point in a list 
                }
                presentCell = cellAnalyzed;
                //Temporal data
                summarizeTemporalAspect(listTimesInCell);

            }
        }
    } // end identifyTimesInCell method

//OK
    /**
     * For each valid temporal interval in the cell, a representative point is
     * defined -- For each valid STI it is created a new representative point,
     * but the computation of -- data summarization (spatial and semantic) it is
     * computed further (on computeCentroid method) -- i.e.
     * summarizeTemporalAspect summarize the Temporal dimension in each cell --
     * then find and define Representative points for cell
     *
     * This method computes the significant temporal intervals (STI) and
     * identifies representative points based on those intervals. These
     * representative points will be used for further computations in the
     * trajectory summarization.
     *
     * @param timeInPoints -- List of temporal information of all points in the
     * cell
     */
    public void summarizeTemporalAspect(ArrayList<Date> timeInPoints) {
        List<STI> significantTemporalIntervals = new ArrayList<>();
        Collections.sort(timeInPoints);

        List<Integer> differences = computeTimeDifferences(timeInPoints);
        float averageDifference = calculateAverage(differences);

        int thresholdDifferences = THRESHOLD_TIME; // Default threshold for <= 2 occurrences
        if (differences.size() > 2) {
            float median = calculateMedian(differences);
            float standardDeviation = calculateStandardDeviation(differences, averageDifference);

            float lowerValue = median - standardDeviation;
            float upperValue = median + standardDeviation;

            List<Integer> validDifferences = removeOutliers(differences, lowerValue, upperValue);
            thresholdDifferences = (int) Math.floor(calculateAverage(validDifferences));
        }

        createTemporalIntervals(timeInPoints, significantTemporalIntervals, thresholdDifferences);

        asortTemporalIntervals(significantTemporalIntervals);

        createRepresentativePoints(significantTemporalIntervals);

    }//end summarizeTemporalAspect method

    /**
     * Compute the representative point of each cell in the spatial grid, last
     * stage of summarization -- all aspects are summarized on it.
     */
    public void computeCentroid() {

        // Ordernate temporal ranking 
        listRepPoint = listRepPoint.stream().sorted().collect(Collectors.toList());

        for (Centroid representativePoint : listRepPoint) {
            resetValuesToSummarization();
            representativeTrajectory.addPoint(representativePoint);
            representativeTrajectory.incrementValue(representativePoint.getPointListSource().size());

            for (Point p : representativePoint.getPointListSource()) {
                // Spatial data
                avgX += p.getX();
                avgY += p.getY();

                //Semantic Data
                Double val;
                String attrActual; //storage attribute order to be used to get the object

                Double[] valuesNumInvalid = {-999.0, -1.0}; //Null values for numerical values

                for (AttributeValue atv : p.getListAttrValues()) {
                    attrActual = "" + atv.getAttibute().getOrder();

//                    System.out.println("att == "+atv.getAttibute());
//                    System.out.println("val == "+atv.getValue());
                    // numeric values - median computation 
                    //in this scope just create bitset with sum and count of values foreach quantitative attribute
                    try {
                        //Vanessa: Algo de errado ao trabalhar com dados numericos

                        val = Double.parseDouble((String) atv.getValue()); // val -1 refers to empty value

                        attributes.get(attributes.indexOf(atv.getAttibute())).setType(SemanticType.NUMERICAL);

                        if (!sematicNumericFusionVal.containsKey(attrActual)) {
                            sematicNumericFusionVal.put(attrActual, new ArrayList<Double>());
                        }
                        //VAnessa: voltar aqui

                        // Add to this key the attribute value if this value is not invalid
                        //i.e., if "val" not is contained in "valuesNumInvalid" 
                        if (!Arrays.asList(valuesNumInvalid).contains(val)) {
                            sematicNumericFusionVal.get(attrActual).add(val);
                        }

                    } catch (java.lang.NumberFormatException e) { //categorical values

                        /*
                            in this scope create the sematicCategoricalSummarizationVal with all possible values of each categorical attribute
                             and add its ids for after this step can computation the frequency of each one,
                             and identify the value more frequency of each qualitative attribute 
                         */
                        if (atv.getAttibute().getType() == null
                                || !atv.getAttibute().getType().equals(SemanticType.NUMERICAL)) {

                            //IF not contains this key - attribute order
                            //Vanessa: update -- key == att
                            attributes.get(attributes.indexOf(atv.getAttibute())).setType(SemanticType.CATEGORICAL);

                            if (!sematicCategoricalSummarizationVal.containsKey(atv.getAttibute())) {
                                sematicCategoricalSummarizationVal.put(atv.getAttibute(), new HashMap<String, Integer>());
                            }

                            // IF this attribute not contains this value
                            if (!sematicCategoricalSummarizationVal.get(atv.getAttibute()).containsKey(atv.getValue())) {
                                //add this value to attribute and initialize the count
                                sematicCategoricalSummarizationVal.get(atv.getAttibute()).put((String) atv.getValue(), 1);
                            } else {
//                            System.out.println("att before: "+sematicCategoricalSummarizationVal.get(atv.getAttibute()));
                                sematicCategoricalSummarizationVal.get(atv.getAttibute()).replace((String) atv.getValue(), sematicCategoricalSummarizationVal.get(atv.getAttibute()).get(atv.getValue()) + 1);
//                            System.out.println("att after: "+sematicCategoricalSummarizationVal.get(atv.getAttibute()));
                            }
                        } // end IF not is numerical

                    }
                } //end FOR of all semantic attributes

            }// end FOR all points for interval time previous analysed 

            // Spatial data
            representativePoint.setSpatialDimension(avgX / representativePoint.getPointListSource().size(), avgY / representativePoint.getPointListSource().size());

            // ---- Semantic data
            //Loop for numeric attributes
            sematicNumericFusionVal.entrySet().forEach((entrada) -> {
                Double median = -999.0;
                Map<Object, Double> newMap = new HashMap<>();

                if (considerNulls) {
//                    
                    //Vanessa: esse IF pode ser retirado
                    if (representativePoint.getPointListSource().size() - entrada.getValue().size()
                            == entrada.getValue().size()) {
                        //When the size of null values and valid values are equals
                        newMap.put("-999.0", 0.5);
                        Collections.sort(entrada.getValue());

                        //Calculates the median value for all numeric attributes of the representative point
                        if (entrada.getValue().size() % 2 == 0) {
                            median = (entrada.getValue().get(entrada.getValue().size() / 2) + entrada.getValue().get(entrada.getValue().size() / 2 - 1)) / 2;
                        } else {
                            median = entrada.getValue().get(entrada.getValue().size() / 2);
                        }
                        newMap.put("" + median, 0.5);

                    } else if ((representativePoint.getPointListSource().size() - entrada.getValue().size())
                            < entrada.getValue().size()) {
                        //When the size of VALID values is more than the size of null values, 
                        //or equal proportion and is not consider null values
                        Collections.sort(entrada.getValue());
                        //Calculates the median value for all numeric attributes of the representative point
                        if (entrada.getValue().size() % 2 == 0) {
                            median = (entrada.getValue().get(entrada.getValue().size() / 2) + entrada.getValue().get(entrada.getValue().size() / 2 - 1)) / 2;
                        } else {
                            median = entrada.getValue().get(entrada.getValue().size() / 2);
                        }
                    }

                } else { // When not consider null values in computation (!considerNulls)
                    //Calculates the median value for all numeric attributes of the representative point 
                    // when the number of valid values is more (or equals) than null values
                    if (!entrada.getValue().isEmpty()) {
                        //Calculates the median value for all numeric attributes of the representative point
                        if (entrada.getValue().size() % 2 == 0) {
                            median = (entrada.getValue().get(entrada.getValue().size() / 2) + entrada.getValue().get(entrada.getValue().size() / 2 - 1)) / 2;
                        } else {
                            median = entrada.getValue().get(entrada.getValue().size() / 2);
                        }

                    }

                }

                //Set attribute value
                if (newMap.isEmpty()) {
                    representativePoint.addAttrValue("" + median,
                            findAttributeForOrder(Integer.parseInt(entrada.getKey())));
                } else {
                    representativePoint.addAttrValue(newMap,
                            findAttributeForOrder(Integer.parseInt(entrada.getKey())));
                }
                //Case of null median not is added (when not consider null values)

            });

            //begin -------- Loop for a categorical attributes
            for (Map.Entry<Object, Map<String, Integer>> allCategorical : sematicCategoricalSummarizationVal.entrySet()) {
                Map<String, Integer> internalCategoricalList
                        = allCategorical.getValue().entrySet()
                                .stream()
                                .sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
                                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                                        (oldValue, newValue) -> oldValue, LinkedHashMap::new));

                System.out.println("Mapa com todos os dados categ: " + allCategorical.getValue().entrySet());
                System.out.println("Pontos do PR: " + representativePoint.getPointListSource().size());
                //Add mode value (tendency) of attribute to representative point

                if (considerNulls) {
                    representativePoint.addAttrValue(normalizeRankingValues(internalCategoricalList,
                            representativePoint.getPointListSource().size(),
                            's'),
                            //                            findAttributeForOrder(Integer.parseInt(allCategorical.getKey()))
                            (SemanticAspect) allCategorical.getKey()
                    );
                } else {
                    representativePoint.addAttrValue(normalizeRankingValuesNotNulls(internalCategoricalList, representativePoint.getPointListSource().size(), 's'),
                            //                            findAttributeForOrder(Integer.parseInt(allCategorical.getKey()))
                            (SemanticAspect) allCategorical.getKey()
                    );
                }
            } // end ------------ Loop for a categorical attributes

            //Reset values to representative point computation
            resetValuesToSummarization();
        }

    }

// OK
    /**
     * Compute the average of the minimum spatial distance of the input MATs
     * points to provide dynamic space segmentation for clustering these input
     * points. Given set of input MATs, with n points, we compute the Euclidean
     * distance d() for each point pi ∈ T with the nearest point pk ∈ T.
     *
     * the spatialThreshold is computed as the average of the minimum spatial --
     * distance the maximun Z value (diagonal size of each cell in the grid) is
     * computed as the max d() of the more distance point
     */
    public void computeSpatialThresholdOutliers() {
        float maxDistanceToZero = 0;
        float auxValueZ;

        float minDistance = Float.MAX_VALUE;
        float sumDistance = 0;
        ArrayList<Float> validDistances = new ArrayList<>();

        for (Point p : points) {
            auxValueZ = Util.euclideanDistanceToZero(p);
            if (auxValueZ > maxDistanceToZero) {
                maxDistanceToZero = auxValueZ;
            }

            for (Point q : points) {
                if (!p.equals(q)) {
                    float localDistance = (float) Util.euclideanDistance(p, q);
                    if (localDistance < minDistance) {
                        minDistance = localDistance;
                    }
                }
            }

            validDistances.add(minDistance);
            sumDistance += minDistance;
            minDistance = Float.MAX_VALUE;
        }

        float avgMinDist = sumDistance / points.size();

        // Remove outliers of minimum spatial distance
        if (validDistances.size() > 1) {
            // Order valid distances
            Collections.sort(validDistances);

            // Compute the median value of the distances
            float medianMinDist = calculateMedian(validDistances);
            // Compute the standard deviation
            float sdMinDist = calculateStandardDeviation(validDistances, avgMinDist);

            // Compute the valid interval
            float lessValueMinDist = medianMinDist - 4 * sdMinDist;
            float upperValueMinDist = medianMinDist + 4 * sdMinDist;

            // Create a new list with valid distances
            ArrayList<Float> validDistancesWithoutOutliers = new ArrayList<>();
            for (float distance : validDistances) {
                if (distance >= lessValueMinDist && distance <= upperValueMinDist && distance != 0.0) {
                    validDistancesWithoutOutliers.add(distance);
                }
            }

            // Compute the spatial threshold -- average of valid values
            spatialThreshold = calculateAverage(validDistancesWithoutOutliers);

            auxMaxZ = maxDistanceToZero / spatialThreshold;
        }
    }

    /**
     * Method core to perform all methods in order to summarize input MATs into
     * one representative MAT.
     *
     * @param file name of file
     * @param ext extension of file
     * @throws IOException
     *
     */
    public void execute(String dir, String file, String ext, String[] lstCategoricalPD, String SEPARATOR, String[] valuesNULL, String[] ignoreColumns, String patternDateInput, float rc, float threshold_rv) throws IOException, ParseException, CloneNotSupportedException {
        initialTemp = new Date();
        //initialization of attribute values (Global attributes according to local data)
        directory = dir;
        filename = file;
        extension = ext;
        this.SEPARATOR = SEPARATOR;
        this.valuesNulls = valuesNULL;
        //Parameter for defining representativeness values and compute relevant cell
        this.threshold_rv = threshold_rv;

        lstCategoricalsPD = Arrays.asList(lstCategoricalPD);
        if (ignoreColumns != null) {
            lstIgnoreCols = Arrays.asList(ignoreColumns);
        }

        //initialization of object of MAT as representative MAT
        representativeTrajectory = new MultipleAspectTrajectory("representative");
        if (!patternDateInput.equals("?")) {
            this.formatDate = new SimpleDateFormat(patternDateInput);
        } else {
            this.representativeTrajectory.setDailyInfo(true);
            dailyInfo = true;
        }

        //initialization of aux attributes
        rId = 0;
        auxTid = "-1";
//        cId = -1;

        //initialization of aux lists
        listTimesInCell = new ArrayList<Date>();
        spatialCellGrid = new HashMap<String, BitSet>();
        sematicNumericFusionVal = new HashMap<String, List<Double>>();
        sematicCategoricalSummarizationVal = new HashMap<Object, Map<String, Integer>>();
        points = new ArrayList<Point>();
        attributes = new ArrayList<SemanticAspect>();
        listTrajectories = new ArrayList<MultipleAspectTrajectory>();
        pointsInCell = new ArrayList<>();

        //aux representative MAT for ordenate
        listRepPoint = new ArrayList<>();

        // Load dataset follow data model representation
        load();

        //rc is defined as the minimun number of points ( calculated by the % of all points) that should have in each cell
        this.rc = rc;
        threshold_rc = rc > 0.0 ? (rc * points.size()) : 2; //If rc is greater than zero sets threshold according with number of points, else sets to 2

        //As auxiliar var to load all input trajectories & all points of trajectories (validation)
        listAllTrajectories = List.copyOf(listTrajectories);
        listTrajectories = new ArrayList<>();
        List<Point> auxClusterPoints = List.copyOf(points);
        points = new ArrayList<Point>();

        loadAllDataset();

        List<MultipleAspectTrajectory> auxDataset = List.copyOf(listTrajectories);
        listTrajectories = null;
        listTrajectories = List.copyOf(listAllTrajectories);
        listAllTrajectories = null;
        listAllTrajectories = List.copyOf(auxDataset);
        points = List.copyOf(auxClusterPoints);

        //######## automation - definition of better Z value - the spatial threshold
        // 1st - Calculates the spatial threshold according with the Z value and point dispersion
        computeSpatialThresholdOutliers();

        // 2nd - Summarize Trajectories into a single representative data
        summarizeTrajetories();

    }

    public void summarizeTrajetories() throws ParseException, CloneNotSupportedException {

        int tempMaxZ = (int) auxMaxZ;
        int tempBetterZ = -1;
        float tempBetterRM = 0;
        float tempZvalueRM = 0;
        float iCoverZ = -1.0f;
        float tempOnlyRM;

        String infoBetterRT = "";
        int count = 0;

        while (tempMaxZ > 1) {

            resetValuesRT();

            cellSizeSpace = (spatialThreshold * tempMaxZ) * 0.7071; // Calcultes size of the cells

            allocateAllPointsInCellSpace(); // Distributes all points in the spatial grid

            identifyTimesInCell(); //Creates the representative trajectory

            computeCentroid();

            if (!representativeTrajectory.getPointList().isEmpty()) {
                tempZvalueRM = (float) medianMeasureRecall();
                System.out.println("Similarity cmputed: " + tempZvalueRM);
                tempOnlyRM = tempZvalueRM;

                iCoverZ = (float) representativeTrajectory.getCoverPoints() / points.size();

                tempZvalueRM = (tempZvalueRM * 0.5f) + (iCoverZ * 0.5f);
//                tempZvalueRM = (tempZvalueRM + iCoverZ) / 2; //*** melhor resultado

                if ((tempZvalueRM * 1.1) >= tempBetterRM) { //***melhor resultado
                    tempBetterZ = tempMaxZ;
                    tempBetterRM = tempZvalueRM;
                    count = 0;
                    betterRT = null;
                    betterRT = (MultipleAspectTrajectory) representativeTrajectory.clone();
                    infoBetterRT = createInfoBetterRT(tempBetterZ, tempBetterRM, tempOnlyRM);

                } else {
                    count++;
                }

            }

            tempMaxZ *= 0.85;

            if (count > 1) {
                break;

            }
        } // fim do laço infinito - // Fim automação
        if (tempBetterZ > 1) {
            for (Map.Entry<String, BitSet> entry : spatialCellGrid.entrySet()) {
                infoBetterRT += entry.getKey()+" | ";
            }
            String outputFile = directory + "output\\" + filename + " rc " + (int) (rc * 100) + " rv " + (int) (threshold_rv * 100) + " - z" + tempBetterZ;
            writeRepresentativeTrajectory(outputFile, infoBetterRT);
            rankMUITAS(outputFile);
            
        }

    }

    public void summarizeTrajetories2() throws ParseException, CloneNotSupportedException {
        int tempBetterZ = -1;
        float tempBetterRM = 0;
        float tempOnlyRM;

        String infoBetterRT = "";
        int count = 0;

        while (auxMaxZ > 1) {
            resetValuesRT();

            cellSizeSpace = (spatialThreshold * auxMaxZ) * 0.7071; // Calculate size of the cells

            allocateAllPointsInCellSpace(); // Distribute all points in the spatial grid

            identifyTimesInCell(); // Create the representative trajectory

            computeCentroid();

            if (!representativeTrajectory.getPointList().isEmpty()) {
                float tempZvalueRM = (float) medianMeasureRecall();
                tempOnlyRM = tempZvalueRM;

                float iCoverZ = (float) representativeTrajectory.getCoverPoints() / points.size();

                tempZvalueRM = (tempZvalueRM + iCoverZ) / 2;

                if ((tempZvalueRM * 1.1) >= tempBetterRM) {
                    tempBetterZ = (int) Math.floor(auxMaxZ);
                    tempBetterRM = tempZvalueRM;
                    count = 0;
                    betterRT = (MultipleAspectTrajectory) representativeTrajectory.clone();
                    infoBetterRT = createInfoBetterRT(tempBetterZ, tempBetterRM, tempOnlyRM);
                } else {
                    count++;
                }
            }

            auxMaxZ *= 0.85;

            if (count > 1) {
                break;
            }
        }
        if (tempBetterZ > 1) {
            String outputFile = directory + "output\\" + filename + " rc " + (int) (rc * 100) + " rv " + (int) (threshold_rv * 100) + " - z" + tempBetterZ;
            writeRepresentativeTrajectory(outputFile, infoBetterRT);
            rankMUITAS(outputFile);
        }

    }

    // #################### 
    // Normalizing our data by ratio value to generate data ranking
    /**
     * For updating the number of occurrences of each rank value by the ratio
     * value. Normalizing the Rank Value Map in the semantic or temporal
     * dimension, where the quantity of occurrences of each attribute value is
     * changed by the ratio of this value in relation to the size of the cell.
     *
     * In temporal dimension the minutes values are converted to valid time
     * information.
     *
     * @param mapRank -- currently the Map of ranking values with number of
     * occurrences for each value
     * @param dimension -- t: temporal and s: semantic
     * @return normalized -- the Map update with ratio values of occurrences.
     */
    public Map<Object, Double> normalizeRankingValues(Map<String, Integer> mapRank, int sizeRP, char dimension) {
        if (mapRank == null || mapRank.isEmpty()) {
            throw new IllegalArgumentException("Invalid input: mapRank must not be null or empty.");
        }
        if (sizeRP <= 0) {
            throw new IllegalArgumentException("Invalid input: sizeRP must be a positive value.");
        }
        int mappedPoints = 0;
        for (Map.Entry<String, Integer> eachValue : mapRank.entrySet()) {
            double trendEachVal = (double) eachValue.getValue() / sizeRP;
            if (trendEachVal >= threshold_rv) {
                mappedPoints += eachValue.getValue();
            } else {
                mapRank.replace(eachValue.getKey(), -1);
            }
        }

        Map<Object, Double> newMap = new HashMap<>();
        for (Map.Entry<String, Integer> eachValue : mapRank.entrySet()) {
            if (eachValue.getValue() != -1) {

                newMap.put(eachValue.getKey(), (double) eachValue.getValue() / mappedPoints);

            }

        }

        Map<Object, Double> newMapSorted = newMap.entrySet().stream()
                .sorted(Map.Entry.<Object, Double>comparingByValue().reversed())
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                        (oldValue, newValue) -> oldValue, LinkedHashMap::new));

        if (dimension == 't') { // temporal dimension
            /*
            In the temporal dimension, the minutes values are converted to valid time according to the predefined mask.
             */
            Map<Object, Double> newTimeMap = new HashMap<>();
            for (Map.Entry<Object, Double> eachInt : newMapSorted.entrySet()) {
                String interval = (String) eachInt.getKey();
                String auxInterval;
                if (interval.contains("-")) {
                    //when temporal dimens. refers to interval
                    auxInterval = formatDateOut.format(Util.convertMinutesToDate(Integer.parseInt(interval.substring(0, interval.indexOf("-")))));
                    auxInterval += " - " + formatDateOut.format(Util.convertMinutesToDate(Integer.parseInt(interval.substring(interval.indexOf("-") + 1))));
                } else {
                    //when it refers to a single occurrence
                    auxInterval = formatDateOut.format(Util.convertMinutesToDate(Integer.parseInt(interval)));
                }

                newTimeMap.put(auxInterval, newMap.get(interval));
            }
            return newTimeMap;
        } else { // Semantic dimension
            return newMapSorted;
        }

    } // end normalizeRankingValues method

    // 
    /**
     * For updating the number of occurrences of each rank value by the ratio
     * value, when user do not want consider null values on generation of
     * representative trajectory. Normalizing the Rank Value Map in the semantic
     * or temporal dimension, where the quantity of occurrences of each
     * attribute value is changed by the ratio of this value in relation to the
     * size of the cell.
     *
     * In temporal dimension the minutes values are converted to valid time
     * information.
     *
     * @param mapRank -- currently the Map of ranking values with number of
     * occurrences for each value
     * @param mappedPoints -- size of origins points of the representativePoint
     * @param dimension -- t: temporal and s: semantic
     * @return normalized -- the Map update with ratio values of occurrences.
     */
    public Map<Object, Double> normalizeRankingValuesNotNulls(Map<String, Integer> mapRank, int mappedPoints, char dimension) {
        if (mapRank == null || mapRank.isEmpty()) {
            throw new IllegalArgumentException("Invalid input: mapRank must not be null or empty.");
        }
        if (mappedPoints <= 0) {
            throw new IllegalArgumentException("Invalid input: mappedPoints must be a positive value.");
        }
        Map<Object, Double> newMap = new HashMap<>();
        double trendEachVal;
        int sizeNotNull = mappedPoints;

        if (mapRank.containsKey(nullValue)) {
            sizeNotNull -= mapRank.get(nullValue);
        }
        if (sizeNotNull < mappedPoints) {
            for (Map.Entry<String, Integer> eachValue : mapRank.entrySet()) {
                if (!eachValue.getKey().equalsIgnoreCase("Unknown")) {
                    trendEachVal = (double) eachValue.getValue() / sizeNotNull;
                    if (trendEachVal >= threshold_rv) {
                        newMap.put(eachValue.getKey(), trendEachVal);
                    }
                }
            }
        }

        Map<Object, Double> newMapSorted = newMap.entrySet().stream()
                .sorted(Map.Entry.<Object, Double>comparingByValue().reversed())
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                        (oldValue, newValue) -> oldValue, LinkedHashMap::new));

        if (dimension == 't') { // temporal dimension
            /*
            In the temporal dimension, the minutes values are converted to valid time according to the predefined mask.
             */
            Map<Object, Double> newTimeMap = new HashMap<>();
            for (Map.Entry<Object, Double> eachInt : newMapSorted.entrySet()) {
                String interval = (String) eachInt.getKey();
                String auxInterval;
                if (interval.contains("-")) {
                    auxInterval = formatDate.format(Util.convertMinutesToDate(Integer.parseInt(interval.substring(0, interval.indexOf("-")))));
                    auxInterval += " - " + formatDate.format(Util.convertMinutesToDate(Integer.parseInt(interval.substring(interval.indexOf("-") + 1))));
                } else {
                    auxInterval = formatDate.format(Util.convertMinutesToDate(Integer.parseInt(interval)));
                }

                newTimeMap.put(auxInterval, newMap.get(interval));
            }
            return newTimeMap;
        } else { // Semantic dimension
            return newMapSorted;
        }

    }

    /**
     * allocate all points of input dataset in spatial cell grid
     */
    public void allocateAllPointsInCellSpace() {
        for (Point p : points) {
            allocateInSpaceCell(p);
        }
    }

    // --- MATSG - auxiliar methods -- to data reset, findAttribute...
    /**
     * Reset all values of all attributes in MAT
     */
    public void resetValuesToSummarization() {
        //Data reset

        //spatial data
        avgX = 0;
        avgY = 0;
//
//        // semantic data (multiple aspects)
        sematicNumericFusionVal.clear();
        sematicCategoricalSummarizationVal.clear();
        presentCell = "";
//
//        //temporal data
        listTimesInCell.clear();
        //all valid points in cell (used for temporal analysis)
        pointsInCell.clear();
//        
    } //end resetValuesToSummarization method

    /**
     * Reset data of values to compute the better RT
     */
    public void resetValuesRT() {
        spatialCellGrid.clear();
        pointsInCell.clear();
        presentCell = null;
        listTimesInCell.clear();
        listRepPoint.clear();
        representativeTrajectory = null;
        representativeTrajectory = new MultipleAspectTrajectory("representative");
        if (dailyInfo == true) {
            representativeTrajectory.setDailyInfo(true);
        }

    }

    /**
     * find the SemanticAspect object by the order
     *
     * @param order
     * @return SemanticAspect
     */
    public SemanticAspect findAttributeForOrder(int order) {
        for (SemanticAspect attr : attributes) {
            if (attr.getOrder() == order) {
                return attr;
            }
        }
        return null;
    } //end findAttributeForOrder method

    // ############## For automatization
    // ################# -- dinamicaly identify better cell size -- #################
    /**
     * Compute the Representativeness Measure by Median for Recall
     *
     * @return Representativeness Measure by Median for Recall
     * @throws ParseException
     */
    public double medianMeasureRecall() throws ParseException {

        if (representativeTrajectory.getPointList().isEmpty()) {
            System.out.println("RT is empty");
            return -1;
        }

//        SimilarityMeasure measure = new SimilarityMeasure();
        MUITAS measure = new MUITAS();
        //Computing and setting thresholds
        //3D with equal weight (0.33) 
        measure.setWeight("SPATIAL", 0.34f);
        measure.setWeight("TIME", 0.33f);

        float auxWeight = 0.33f / (attributes.size());

        //System.out.println("Atributes: " + attributes);
        for (SemanticAspect eachAtt : attributes) {
            measure.setWeight(eachAtt, auxWeight);
            if (eachAtt.getType().equals(SemanticType.NUMERICAL)) {
                measure.setThreshold(eachAtt, 10);
            }
        }

        //add spatial threshold -- Test: spatialThreshold x 2 -- pensando em atingir a distância de até 2 células
        measure.setThreshold("SPATIAL", (spatialThreshold * 2));

        double repMeasure = 0;
        List<Double> listValues = new ArrayList<>();

        for (MultipleAspectTrajectory eachTraj : listTrajectories) {
            System.out.println("###");
            System.out.println("Computing similarity");

            listValues.add(measure.similarityOf(representativeTrajectory, eachTraj));
            System.out.println("Similarity: " + listValues.get(listValues.size() - 1));

        }
        //after computed measure with each T and RT, it is computed median value
        repMeasure = calculateMedian(listValues);

        return repMeasure;

    }

    // ############## For validation #############
    private void loadAllDataset() throws IOException, ParseException {

        java.io.Reader input = new FileReader(directory + filenameFullDataset + extension);
        BufferedReader reader = new BufferedReader(input);

        String datasetRow = reader.readLine();

        datasetRow = reader.readLine();
        String[] datasetColumns;
        //EoF - To get the trajectory data of dataset of each line
        while (datasetRow != null) {
            datasetColumns = datasetRow.toUpperCase().split(SEPARATOR);
            addAttributeValues(datasetColumns);
            datasetRow = reader.readLine();
        }

        reader.close();

    }

    public void rankMUITAS(String fileOutput) throws ParseException {
        if (betterRT.getPointList().isEmpty()) {
            System.out.println("RT is empty");
            return;
        }

        MUITAS measure = new MUITAS();

        // Compute weights
        float totalWeight = 1.0f;
        float spatialWeight = 0.34f;
        float timeWeight = 0.33f;
        float attributeWeight = totalWeight / attributes.size();
        measure.setWeight("SPATIAL", spatialWeight);
        measure.setWeight("TIME", timeWeight);
        for (SemanticAspect eachAtt : attributes) {
            measure.setWeight(eachAtt, attributeWeight);
            if (eachAtt.getType().equals(SemanticType.NUMERICAL)) {
                measure.setThreshold(eachAtt, 10);
            }
        }

        // Set spatial threshold
        measure.setThreshold("SPATIAL", spatialThreshold * 2);

        // Compute rank measures
        Map<MultipleAspectTrajectory, Double> trajectoryRanks = new HashMap<>();
        StringBuilder infoMeasure = new StringBuilder();
        int countPrecisionRetrieved = 0;
        int countRecallRetrieved = 0;

        for (MultipleAspectTrajectory eachTraj : listAllTrajectories) {
            double similarity = measure.similarityOf(betterRT, eachTraj);
            trajectoryRanks.put(eachTraj, similarity);
        }

        trajectoryRanks = trajectoryRanks.entrySet().stream()
                .sorted(Map.Entry.comparingByValue(Comparator.reverseOrder()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                        (oldValue, newValue) -> oldValue, LinkedHashMap::new));

        int countTClass = 0;
        for (Map.Entry<MultipleAspectTrajectory, Double> retrievedT : trajectoryRanks.entrySet()) {
            countPrecisionRetrieved++;
            MultipleAspectTrajectory retrievedTrajectory = retrievedT.getKey();
            infoMeasure.append(retrievedTrajectory.getId()).append(", ")
                    .append(retrievedT.getValue()).append(", ")
                    .append(countPrecisionRetrieved).append(", ");

            if (listTrajectories.contains(retrievedTrajectory)) {
                countTClass++;
                infoMeasure.append("1");

                if (countPrecisionRetrieved <= listTrajectories.size()) {
                    countRecallRetrieved++;
                }
            } else {
                infoMeasure.append("0");
            }
            infoMeasure.append("\n");

            if (countTClass == listTrajectories.size()) {
                break;
            }
        }

        // Write validation results to CSV file
        try {
            CSVWriter mxWriter = new CSVWriter(fileOutput + "[Validation]" + extension);
            mxWriter.writeLine("Method validation information:");
            mxWriter.writeLine("|Ground Truth|, |all input dataset|, |T.P.retrieved|, Precision, |T.R.retrieved|, Recall, F-Score");
            mxWriter.writeLine(listTrajectories.size() + ", " + listAllTrajectories.size() + ", "
                    + countPrecisionRetrieved + ", " + formatNumber.format((double) listTrajectories.size() / countPrecisionRetrieved) + ", "
                    + countRecallRetrieved + ", " + formatNumber.format((float) countRecallRetrieved / listTrajectories.size()) + ", ??");
            mxWriter.writeLine("##");
            mxWriter.writeLine("Measure description:");
            mxWriter.writeLine("Trajectory ID, MUITAS, #rank, Ground Truth?");
            mxWriter.writeLine(infoMeasure.toString());
            mxWriter.flush();
            mxWriter.close();

            String fileCompleteValidation = directory + "output\\" + filename + "[Validation]" + extension;
            CSVWriter valWriter;
            if (!new File(fileCompleteValidation).exists()) {
                valWriter = new CSVWriter(fileCompleteValidation);
                valWriter.writeLine("Method validation information:");
                valWriter.writeLine("Setting rv, Setting rc, |Ground Truth|, |all input dataset|, |T.P.retrieved|, Precision, |T.R.retrieved|, Recall, F-Score");
            } else {
                valWriter = new CSVWriter(fileCompleteValidation, true);
            }
            valWriter.writeLine(formatNumber.format(threshold_rv) + ", " + formatNumber.format(threshold_rc / points.size()) + ", "
                    + listTrajectories.size() + ", " + listAllTrajectories.size() + ", "
                    + countPrecisionRetrieved + ", " + formatNumber.format((double) listTrajectories.size() / countPrecisionRetrieved) + ", "
                    + countRecallRetrieved + ", " + formatNumber.format((float) countRecallRetrieved / listTrajectories.size()) + ", ??");
            valWriter.flush();
            valWriter.close();

        } catch (IOException e) {
            throw new RuntimeException("Error occurred while ranking input trajectories against RT", e);
        }
    }

    // Reusable codes
    /**
     * calculates the average of a list of integers correctly accept "float" or
     * "integer" values
     *
     * @param values
     * @return average value
     */
    private <T extends Number> float calculateAverage(List<T> values) {
        float sum = 0;
        for (T val : values) {
            sum += val.floatValue();

        }
        return (float) sum / values.size();
    }

    /**
     * calculates the median of a list of integers correctly.
     *
     * @param values
     * @return median value
     */
    private <T extends Number & Comparable<? super T>> float calculateMedian(List<T> values) {
        Collections.sort(values, Comparator.naturalOrder());
        int size = values.size();
        int middleIndex = size / 2;
        if (size % 2 == 1) {
            return values.get(middleIndex).floatValue();
        } else {
            T value1 = values.get(middleIndex - 1);
            T value2 = values.get(middleIndex);
            return (value1.floatValue() + value2.floatValue()) / 2;
        }
    }

    private <T extends Number> float calculateStandardDeviation(List<T> values, float average) {
        float sumOfSquares = 0;
        for (T value : values) {
            float differenceMinusAverage = value.floatValue() - average;
            sumOfSquares += differenceMinusAverage * differenceMinusAverage;
        }
        float meanOfSquares = sumOfSquares / values.size();
        return (float) Math.sqrt(meanOfSquares);
    }

    private List<Integer> removeOutliers(List<Integer> values, float lowerValue, float upperValue) {
        List<Integer> validValues = new ArrayList<>();
        for (int val : values) {
            if (val >= lowerValue && val <= upperValue) {
                validValues.add(val);
            }
        }
        return validValues;
    }

    private void createTemporalIntervals(List<Date> timeInPoints, List<STI> significantTemporalIntervals, int threshold) {
        int count = 1;
        TemporalAspect currentInterval = null;

        for (int i = 0; i < timeInPoints.size(); i++) {
            if (i != timeInPoints.size() - 1
                    && (TimeUnit.MINUTES.convert(timeInPoints.get(i).getTime(), TimeUnit.MILLISECONDS) + threshold)
                    >= (TimeUnit.MINUTES.convert(timeInPoints.get(i + 1).getTime(), TimeUnit.MILLISECONDS))) {
                if (currentInterval == null) {
                    currentInterval = new TemporalAspect(timeInPoints.get(i));
                }
                count++;
            } else {
                if (currentInterval == null) {
                    currentInterval = new TemporalAspect(timeInPoints.get(i));
                } else {
                    currentInterval.setEndTime(timeInPoints.get(i));
                }
                significantTemporalIntervals.add(new STI(currentInterval, (float) count / timeInPoints.size()));
                count = 1;
                currentInterval = null;
            }
        }
    }

    private void asortTemporalIntervals(List<STI> significantTemporalIntervals) {
        significantTemporalIntervals.sort(Comparator.comparing(STI::getProportion).reversed());
    }

    private void createRepresentativePoints(List<STI> significantTemporalIntervals) {
        for (STI interval : significantTemporalIntervals) {
            if (interval.getProportion() >= threshold_rv) {
                Centroid representativePoint = createCentroidForSTI(interval);

                if (!representativePoint.getPointListSource().isEmpty()) {
                    representativePoint.setSti(interval);
                    listRepPoint.add(representativePoint);
                }
            }
        }
    }

    private Centroid createCentroidForSTI(STI eachSTI) {
        Centroid tempRepP = new Centroid();
        for (Point p : pointsInCell) {
            if (eachSTI.getInterval().isInInterval(p.getTime().getStartTime())) {
                tempRepP.addPoint(p);
            }
        }
        tempRepP.setCellReference(presentCell);
        return tempRepP;
    }

    private List<Integer> computeTimeDifferences(ArrayList<Date> timeInPoints) {
        List<Integer> differences = new ArrayList<>();
        for (int i = 1; i < timeInPoints.size(); i++) {
            int auxDif = (int) TimeUnit.MINUTES.convert(timeInPoints.get(i).getTime() - timeInPoints.get(i - 1).getTime(), TimeUnit.MILLISECONDS);
            if (auxDif > 0) {
                differences.add(auxDif);
            }
        }
        return differences;
    }

// ----- Print data -- typing
    /**
     * Writes the generated representative trajectory in a new .csv file
     *
     * @param fileOutput -- output file name
     */
    public void writeRepresentativeTrajectory(String fileOutput, String infoBetterRT) {
        try {
            CSVWriter mxWriter = new CSVWriter(fileOutput + extension);
            mxWriter.writeLine("Method runtime information:");
            mxWriter.writeLine("Start timestamp: " + initialTemp);
            mxWriter.writeLine("End timestamp: " + new Date());
            mxWriter.writeLine("##");
            mxWriter.writeLine("Info input dataset:");
            mxWriter.writeLine("|input.T|, |input.T.points|");
            mxWriter.writeLine(listTrajectories.size() + ", " + points.size());
            mxWriter.writeLine("##");
            mxWriter.writeLine("RT setting infos:");
            mxWriter.writeLine("thresholdCellSize, |rt|, CellSize, "
                    + "tauRelevantCell, minPointsRC, "
                    + "tauRepresentativenessValue, |cell|, "
                    + "RepresentativenessMeasure, |cover RT|");
            mxWriter.writeLine(infoBetterRT);
            mxWriter.writeLine("##");
            mxWriter.writeLine("RT description:");
            for (Point p : betterRT.getPointList()) {
                mxWriter.writeLine(p.toString());
                mxWriter.flush();
            }
            mxWriter.flush();
            mxWriter.close();
        } catch (IOException e) {
            System.err.println("Error on writing RT: " + e.toString());
        }
    }

    private String createInfoBetterRT(int tempBetterZ, float tempBetterRM, float tempOnlyRM) {
        return tempBetterZ + ", "
                + betterRT.getPointList().size() + ", "
                + cellSizeSpace + ", " + rc + ", "
                + threshold_rc + ", " + threshold_rv + ", "
                + spatialCellGrid.size() + ", "
                + tempBetterRM + ", "
                + betterRT.getCoverPoints() + ", "
                + tempOnlyRM;
    }

    //################# getter and setter
    // Setter methods
    public void setSEPARATOR(String SEPARATOR) {
        this.SEPARATOR = SEPARATOR;
    }

    public void setINDEX_SEMANTIC(int INDEX_SEMANTIC) {
        this.INDEX_SEMANTIC = INDEX_SEMANTIC;
    }

    public void notConsiderNulls() {
        considerNulls = false;
    }

    public void setFilenameFullDataset(String filenameFullDataset) {
        this.filenameFullDataset = filenameFullDataset;
    }

}
