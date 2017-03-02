import java.util.ArrayList;
import java.util.List;

public class Burst_detection {
	double[] q;
	double[] bursts;//突发状态标为1，
	bursts_range br;
	List<bursts_range> lbr;
	double p[] ;
	
	public Burst_detection() {
		lbr=new ArrayList<bursts_range>();
	}
    public static double getBinomial(double N, double k)
    {
        double m = N-k;
        double min = k;
        double max = m;
        double t = 0;

        double NN=1;
        double kk=1;
        if(min>max)
        {
            t=min;
            min = max;
            max=t;
        }
        while(N>max)
        {
            NN=NN*N;
            N--;
        }
        while(min>0)
        {
            kk=kk*min;
            min--;
        }
        return NN/kk;
    }
 
	
	
	public static int where(double[] arr) {

		int i;
		double minValue = arr[0];
		int j = 0;
		for (i = 0; i < arr.length; i++) {
			if (arr[i] < minValue) 
				j = i;
			minValue = arr[i];
		}
		return j;
	}

	public static double nansum(double[] array) {
		double sum = 0;
		for (int j = 0; j < array.length; j++) {
			sum += array[j];
		}

		return sum;
	}

	public double tau(double i1, double i2, double gamma, double n) {
		if (i1 >= i2) {
			return 0;
		} else {
			return (i2 - i1) * gamma * Math.log(n);
		}
	}

	public double fit(double d, double r, double p) {
		return -Math.log((double) getBinomial(d, r) * (Math.pow(p, r))
				* (Math.pow((1 - p), (d - r))));
	}

	public void burst_detection(double[] r, double[] d, int n, double s,
			double gamma, int smooth_win) {
		int k = 2;
		int size_r = r.length;
		double[] temp_p = new double[size_r];
		double[] temp_p_temp = new double[size_r];
		int real_n = 0;
		if (smooth_win > 1) {
			real_n = 0;
			for (int i = 0; i < size_r; i++)
				temp_p_temp[i] = r[i] / d[i];
			for (int i = 0; i < size_r; i++)
				temp_p[i] = -1;
			for (int i = 0; i < (size_r - smooth_win + 1); i++) {
				double temp_sum = 0;
				for (int j = 0; j < smooth_win; j++) {
					temp_sum += temp_p_temp[i];
				}
				temp_p[(i + 1) / 2] = temp_sum / smooth_win;
			}

			for (int i = 0; i < size_r; i++) {
				if ((i > ((smooth_win + 1) / 2))
						&& ((i < (size_r - smooth_win + 1)))) {
					r[i] = (int) (temp_p[i] * d[i]);
					real_n++;
				} else {
					r[i] = -1;
				}
			}

		} else {
			real_n = n;
		}

		p = new double[2];

		p[0] = nansum(r) / nansum(d);
		p[1] = p[0] * s;
		if (p[1] > 1)
			p[1] = 0.99999;

		double[][] cost = new double[n][k];		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < k; j++) {
				cost[i][j]=-1;
			}
		}
		 q= new double[n];

		int t = 0;
		for (t = (smooth_win - 1) / 2; t < (((smooth_win - 1) / 2) + real_n); t++) {
			for (int j = 0; j < k; j++) {
				if (t == ((smooth_win - 1) / 2)) {
					cost[t][j] = fit(d[t], r[t], p[j]);
				} else {
					cost[t][j] = tau(q[t - 1], j, gamma, real_n)
							+ fit(d[t], r[t], p[j]);
				}
			}
			q[t] = where(cost[t]);
		}

		
		for (int i = 0; i < n; i++) {
			System.out.println(q[i]);
		}
	}

	
	public void enumerate_bursts(int n){
		bursts=new double[n];
		for (int i = 0; i < n; i++) {
			bursts[i]=0;
		}
		boolean burst = false;
		for (int i = 1; i < q.length; i++) {
			
			if((burst == false)&(q[i]>q[i-1]))
				{	
					br=new bursts_range();
					br.begin=i;
					burst = true;
				}
			if((burst == true)&(q[i]<q[i-1]))
				{	
					br.end=i;
					lbr.add(br);
					burst = false;
				}
		}
		
	}
	
	
	public void burst_weights(int n,double[] r, double[] d){
		double  cost_diff_sum ;
		for (int i = 0; i < lbr.size(); i++) {
			cost_diff_sum = 0;
			bursts_range br_tmp=lbr.get(i);
			if ((br_tmp.begin!=-1)&&(br_tmp.end!=-1)) {
				for (int j =br_tmp.begin; j < br_tmp.end; j++) {
					cost_diff_sum = cost_diff_sum + (fit(d[j],r[j],p[0]) - fit(d[j],r[j],p[1]));
				}				
			}
		System.out.println(cost_diff_sum);
		}
	}
	
	
	public static void main(String[] args) {

		double[] r = new double[] { 0, 2, 1, 6, 7, 2, 8, 7, 2, 1 };
		double[] d = new double[] { 9, 11, 12, 10, 10, 8, 12, 10, 13, 11 };
		int n = r.length;
		double s = 2;
		double gamma = 1;
		int smooth_win = 1;
		
		Burst_detection bd = new Burst_detection();
		bd.burst_detection(r, d, n, s, gamma, smooth_win);
		bd.enumerate_bursts(n);
		bd.burst_weights(n,r,d);
	}
}
