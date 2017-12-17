public class MatrixCalculation {

    private int formula;
    private int orderOfFunction;
    private double C;
    private double minError;
    private MRRunner runner;

    MatrixCalculation(int formula, int order, double c, double minError) {
        this.formula = formula;
        this.orderOfFunction = order;
        this.C = c;
        this.minError = minError;
        this.runner = new MRRunner();
    }

    public double[][] calculateMatrix(double betta, double gamma) {
        double[][] matrix = new double[orderOfFunction][];
        runner.initJacobiCalc(formula, orderOfFunction);

        for (int k = 1; k <= orderOfFunction; k++) {
            double deltaTao = Math.sqrt(8 * minError / Math.abs(derivativeMax(k, betta, gamma)));

            double first = formula(k, betta, 0, gamma);
            double current = formula(k, betta, deltaTao, gamma);
            double nextTao = 2 * deltaTao;

            int countSmaller = 0;
            int smallerThreshold = 500;
            while (true) {
                double next = formula(k, betta, nextTao, gamma);
                if (Math.signum(current - first) != Math.signum(next - current)) {
                    if (Math.abs(current) < minError) break;
                }
                if (Math.abs(current) < minError) {
                    countSmaller++;
                    if (countSmaller >= smallerThreshold) break;
                } else {
                    countSmaller = 0;
                }

                first = current;
                current = next;
                nextTao += deltaTao;
            }

            /*//todo: remove
            int N = 10;
            deltaTao = 0.05;
            //todo: remove*/
            int N = (int) Math.floor(nextTao / deltaTao + 0.5);

            runner.addInputParamsForKFuncOrder(formula, k, N, gamma, betta, deltaTao);
         /*   matrix[k - 1] = new double[N];
            for (int i = 0; i < N; i++) {
                double tao = deltaTao * i;
                matrix[k - 1][i] = formula(k, betta, tao, gamma);
            }*/
        }

        matrix = runner.calcJacobiFunction();

        return matrix;
    }

    public double[][] transposeMatrix(double[][] matrix) {
        return runner.transposeMatrix(matrix);
    }

    public double[][] multiplyMatrixByTransposed(double[][] matrix) {
        return runner.multiplyMatrixByTransposed(matrix);
    }

    private double derivate(int k, double betta, double tao, double gamma12) {
        double sum = 0;
        double sign = 1;
        for (int s = 0; s <= k; s++) {
            sum += sign * C(k, s) * doubleC(k + s + betta, s) * pow(2 * s + 1, 1) * Math.exp(-(2 * s + 1) * C * gamma12 * tao / 2);

            sign *= -1;
        }
        return pow(-C * gamma12 / 2, 1) * sum;
    }

    private double p(int k, double betta1, double tao, double gamma1) {
        double res = 0;
        int sign = 1;
        for (int s = 0; s <= k; s++) {
            res += C(k, s) * doubleC(k + s + betta1, s) * sign * Math.exp(-(2 * s + 1) * 2 * gamma1 * tao / 2);

            sign *= -1;
        }
        return res;
    }

    private double integral(int k, double betta13, double tao, double gamma13) {
        int n = 1;
        double res = 0;
        double sign = 1;
        for (int s = 0; s <= k; s++) {
            double binoms = C(k, s) * doubleC(k + s + betta13, s);
            double exps = sign * Math.exp(-(2 * s + 1) * C * gamma13 * tao / 2);
            double rowSum = 0;
            for (int j = 0; j <= n; j++) {
                rowSum += f(n) * pow(tao, n - j) / (f(n - j) * pow(C * gamma13 * (2 * s + 1) / 2, j + 1));
            }
            res += binoms * exps * rowSum;
            sign *= -1;
        }
        return res;
    }

    private double formula(int k, double betta, double tao, double gamma) {
        if (formula == 0) return derivate(k, betta, tao, gamma);
        if (formula == 1) return p(k, betta, tao, gamma);
        if (formula == 2) return integral(k, betta, tao, gamma);
        return -1;
    }

    //formula 1.65
    private double derivativeMax(int k, double betta, double gamma) {
        return derivativeCommon(2, k, betta, 0, gamma);
    }

    private double derivativeCommon(int n, int k, double betta, double tao, double gamma) {
        double sum = 0;
        double sign = 1;
        for (int s = 0; s <= k; s++) {
            sum += sign * c(k, s) * c(k + s + betta, s) * power(2 * s + 1, n) * Math.exp(-(2 * s + 1) * C * gamma * tao / 2);
            sign *= -1;
        }
        return power(-C * gamma / 2, n) * sum;
    }

    private double c(double n, double k) {
        return Gamma.gamma(n + 1) / Gamma.gamma(k + 1) / Gamma.gamma(n - k + 1);
    }

    private double power(double x, int n) {
        return Math.pow(x, n);
    }

    private static double[][] cacheC = new double[100][100];

    static double C(int n, int k) {
        if (cacheC[n][k] != 0) return cacheC[n][k];
        return cacheC[n][k] = f(n) / f(k) / f(n - k);
    }

    static double doubleC(double n, double k) {
        return Gamma.gamma(n + 1) / Gamma.gamma(k + 1) / Gamma.gamma(n - k + 1);
        //return Gamma.gamma(n) / Gamma.gamma(k) / Gamma.gamma(n - k);
    }

    static double f(int n) {
        if (n == 0) return 1;
        return f(n - 1) * n;
    }

    static double pow(double x, int n) {
        if (n == 0) return 1;
        if (n == 1) return x;
        if (n == 2) return x * x;
        return Math.pow(x, n);
    }
}