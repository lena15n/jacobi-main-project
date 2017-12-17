import java.io.*;
import java.util.*;

public class MRRunner {
    private String rootDir = System.getProperty("user.dir").concat("/");
    private String inputFile;
    private String outputFile;

    private int formula;
    private int orderOfFunction;

    public void initJacobiCalc(int formula, int orderOfFunction) {
        this.formula = formula;
        this.orderOfFunction = orderOfFunction;
    }

    public void initFiles() {
        if (inputFile == null) {
            inputFile = "input_jacobi";
            File inputF = new File(rootDir + inputFile);
            try {
                inputF.createNewFile();

            } catch (IOException e) {
                System.err.println(e.getMessage());
            }
        }

        outputFile = "out_jacobi";
    }

    public void addInputParamsForKFuncOrder(int formula, int k, int n, double gamma, double betta, double deltaTao) {
        try {
            PrintWriter printWriter = new PrintWriter(new BufferedWriter(
                    new FileWriter(rootDir + inputFile, true)));
            printWriter.println(String.format("\t%d %d %d %.15f %.15f %.15f", formula, k, n, gamma, betta, deltaTao));
            printWriter.close();

        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
    }

    public double[][] calcJacobiFunction() {
        return runMRJob("src/main/resources/jacobi.jar", null,
                new JacobiMatrixBuilder(orderOfFunction));
    }

    public double[][] transposeMatrix(double[][] matrix) {
        initFiles();
        DefaultMatrixBuilder builder = new DefaultMatrixBuilder(matrix[0].length, matrix.length);
        builder.writeMatrix(matrix, inputFile);
        String rowsCountParam = String.valueOf(matrix.length);
        return runMRJob("src/main/resources/transp.jar", new String[]{rowsCountParam}, builder);
    }

    public double[][] multiplyMatrixByTransposed(double[][] matrix) {
        initFiles();
        MultiplyMatrixBuilder builder = new MultiplyMatrixBuilder(matrix.length, matrix.length);
        builder.writeMatrix(matrix, inputFile);
        String rowsCountParam = String.valueOf(matrix.length);
        String columnsCountParam = String.valueOf(matrix[0].length);
        return runMRJob("src/main/resources/transp-mult-invert.jar",
                new String[]{rowsCountParam, columnsCountParam}, builder);
    }

    class MultiplyMatrixBuilder extends DefaultMatrixBuilder {
        MultiplyMatrixBuilder(int rows, int columns) {
            super(rows, columns);
        }

        @Override
        public double[][] buildMatrix() {
            double[][] matrix = new double[outRows][];
            for (int i = 0; i < matrix.length; i++) {
                matrix[i] = new double[outColumns];
            }

            File file = new File(rootDir + outputFile);
            try {
                Scanner sc = new Scanner(new FileInputStream(file));
                while (sc.hasNext()) {
                    String[] line = sc.nextLine().split("\\t");
                    String[] indices = line[0].split(" ");
                    int i = Integer.valueOf(indices[0]);
                    int j = Integer.valueOf(indices[1]);

                    matrix[i][j] = Double.valueOf(line[1]);
                }

            } catch (FileNotFoundException e) {
                System.err.println(e.getMessage());
            }

            return matrix;
        }
    }

    class DefaultMatrixBuilder implements MatrixBuilder {
        protected int outRows;
        protected int outColumns;

        DefaultMatrixBuilder(int outRows, int outColumns) {
            this.outRows = outRows;
            this.outColumns = outColumns;
        }

        @Override
        public double[][] buildMatrix() {
            double[][] matrix = new double[outRows][];
            File file = new File(rootDir + outputFile);
            try {
                Scanner sc = new Scanner(new FileInputStream(file));
                while (sc.hasNext()) {
                    String[] line = sc.nextLine().split("\\t");
                    int i = Integer.valueOf(line[0]);

                    StringTokenizer value = new StringTokenizer(line[1]);

                    matrix[i] = new double[outColumns];
                    for (int j = 0; j < matrix[0].length; j++) {
                        matrix[i][j] = Double.valueOf(value.nextToken());
                    }
                }

            } catch (FileNotFoundException e) {
                System.err.println(e.getMessage());
            }

            return matrix;
        }

        public void writeMatrix(double[][] matrix, String inputFile) {
            File file = new File(inputFile);
            try {
                PrintWriter pw = new PrintWriter(file);
                for (int i = 0; i < matrix.length; i++) {
                    pw.print(i + "\t");
                    for (int j = 0; j < matrix[0].length; j++) {
                        pw.print(String.format("%.15f", matrix[i][j]));
                        pw.print(" ");
                    }
                    pw.println();
                }
                pw.close();

            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }

    class JacobiMatrixBuilder implements MatrixBuilder {
        private int rows;

        JacobiMatrixBuilder(int rows) {
            this.rows = rows;
        }

        public double[][] buildMatrix() {
            double[][] matrix = new double[rows][];
            File file = new File(rootDir + outputFile);
            try {
                Scanner sc = new Scanner(new FileInputStream(file));
                while (sc.hasNext()) {
                    String[] line = sc.nextLine().split("\\t");
                    String[] key = line[0].split(" ");
                    int k = Integer.valueOf(key[0]);
                    int i = Integer.valueOf(key[1]);

                    String[] value = line[1].split(" ");
                    double function = Double.valueOf(value[0]);

                    if (matrix[k - 1] == null) {
                        int n = Integer.valueOf(value[1]);
                        matrix[k - 1] = new double[n];
                    }

                    matrix[k - 1][i] = function;
                }

            } catch (FileNotFoundException e) {
                System.err.println(e.getMessage());
            }

            return matrix;
        }
    }

    private interface MatrixBuilder {
        double[][] buildMatrix();
    }

    private double[][] runMRJob(String pathToJobJar, String[] args, MatrixBuilder builder) {
        try {
            if (inputFile == null) {
                initFiles();
            }

            List<String> command = new ArrayList<>();
            command.add("hadoop");
            command.add("jar");
            command.add(pathToJobJar);
            command.add(inputFile);
            command.add(outputFile);
            if (args != null) {
                command.addAll(Arrays.asList(args));
            }

            ProcessBuilder pb = new ProcessBuilder(command)
                    .redirectErrorStream(true)
                    // Inherit System.out as redirect output stream
                    .redirectOutput(ProcessBuilder.Redirect.INHERIT);

            long time = System.currentTimeMillis();

            Process p = pb.start();
            p.waitFor();

            time = System.currentTimeMillis() - time;
            System.out.println("Time: " + time);

            double[][] result = builder.buildMatrix();

            File inFile = new File(rootDir + inputFile);
            File outFile = new File(rootDir + outputFile);
            inFile.delete();
            outFile.delete();

            return result;

        } catch (IOException | InterruptedException e) {
            System.err.println(e.getMessage());
        }

        return null;
    }
}