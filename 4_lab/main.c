#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

struct particle {
	float x, y, z;
};

const float G = 6.67e-11;

void calculate_forces_var1(struct particle *p, struct particle *f, float *m, int n, int x)
{
	#pragma omp parallel num_threads(x)
	{
		#pragma omp for 
		for (int i = 0; i < n - 1; i++) {
			for (int j = i + 1; j < n; j++) {
				float dist = sqrtf(powf(p[i].x - p[j].x, 2) + 
									powf(p[i].y - p[j].y, 2) + 
									powf(p[i].z - p[j].z, 2));
				float mag = (G * m[i] * m[j]) / powf(dist, 2);
				struct particle dir = {
					.x = p[j].x - p[i].x,
					.y = p[j].y - p[i].y,
					.z = p[j].z - p[i].z
				};

				#pragma omp critical
				{
					f[i].x += mag * dir.x / dist;
					f[i].y += mag * dir.y / dist;
					f[i].z += mag * dir.z / dist;
					f[j].x -= mag * dir.x / dist;
					f[j].y -= mag * dir.y / dist;
					f[j].z -= mag * dir.z / dist;
				}
				
			}
		}
	}
	
}

void calculate_forces_var2(struct particle *p, struct particle *f, float *m, int n)
{
	#pragma omp for schedule(dynamic, 4) nowait
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			float dist = sqrtf(powf(p[i].x - p[j].x, 2) + 
								powf(p[i].y - p[j].y, 2) + 
								powf(p[i].z - p[j].z, 2));
			float mag = (G * m[i] * m[j]) / powf(dist, 2);
			struct particle dir = {
				.x = p[j].x - p[i].x,
				.y = p[j].y - p[i].y,
				.z = p[j].z - p[i].z
			};

			#pragma omp atomic
			f[i].x += mag * dir.x / dist;
			#pragma omp atomic
			f[i].y += mag * dir.y / dist;
			#pragma omp atomic
			f[i].z += mag * dir.z / dist;
			#pragma omp atomic
			f[j].x -= mag * dir.x / dist;
			#pragma omp atomic
			f[j].y -= mag * dir.y / dist;
			#pragma omp atomic
			f[j].z -= mag * dir.z / dist;
		}
	}
}

omp_lock_t *locks; 

void calculate_forces_var3(struct particle *p, struct particle *f, float *m, int n)
{
	#pragma omp for schedule(dynamic, 4) nowait
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			float dist = sqrtf(powf(p[i].x - p[j].x, 2) + 
								powf(p[i].y - p[j].y, 2) + 
								powf(p[i].z - p[j].z, 2));
			float mag = (G * m[i] * m[j]) / powf(dist, 2);
			struct particle dir = {
				.x = p[j].x - p[i].x,
				.y = p[j].y - p[i].y,
				.z = p[j].z - p[i].z
			};

			omp_set_lock(&locks[i]);
			f[i].x += mag * dir.x / dist;
			f[i].y += mag * dir.y / dist;
			f[i].z += mag * dir.z / dist;
			omp_unset_lock(&locks[i]);
			omp_set_lock(&locks[j]);
			f[j].x -= mag * dir.x / dist;
			f[j].y -= mag * dir.y / dist;
			f[j].z -= mag * dir.z / dist;
			omp_unset_lock(&locks[j]);
		}
	}
}

void calculate_forces_var4(struct particle *p, struct particle *f, float *m, int n)
{
	#pragma omp for schedule(dynamic, 4) nowait
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				continue;
			}
			float dist = sqrtf(powf(p[i].x - p[j].x, 2) + 
								powf(p[i].y - p[j].y, 2) + 
								powf(p[i].z - p[j].z, 2));
			float mag = (G * m[i] * m[j]) / powf(dist, 2);
			struct particle dir = {
				.x = p[j].x - p[i].x,
				.y = p[j].y - p[i].y,
				.z = p[j].z - p[i].z
			};

			f[i].x += mag * dir.x / dist;
			f[i].y += mag * dir.y / dist;
			f[i].z += mag * dir.z / dist;
		}
	}
}

void calculate_forces_var5(struct particle *p, struct particle *f[], float *m, int n)
{
	int tid = omp_get_thread_num();
	int nthreads = omp_get_num_threads();
	for (int i = 0; i < n; i++) {
		f[tid][i].x = 0;
		f[tid][i].y = 0;
		f[tid][i].z = 0;
	}
	// #pragma omp for collapse(2) schedule(dynamic, 1024) nowait
	// for (int i = 0; i < n - 1; i++) {
	// 	for (int j = i + 1; j < n; j++) {
	#pragma omp for schedule(dynamic, 25) nowait
	for (int w = 0; w < n * (n - 1) / 2; w++) {
		int i = (-1 + sqrt(1.0 + 8.0 * w)) / 2;
		int j = w - i * (i + 1) / 2;

		float dist = sqrtf(powf(p[i].x - p[j].x, 2) + 
							powf(p[i].y - p[j].y, 2) + 
							powf(p[i].z - p[j].z, 2));

		float mag = (G * m[i] * m[j]) / powf(dist, 2);

		struct particle dir = { 
			.x = p[j].x - p[i].x, 
			.y = p[j].y - p[i].y, 
			.z = p[j].z - p[i].z 
		};

		f[tid][i].x += mag * dir.x / dist;
		f[tid][i].y += mag * dir.y / dist;
		f[tid][i].z += mag * dir.z / dist;
		f[tid][j].x -= mag * dir.x / dist;
		f[tid][j].y -= mag * dir.y / dist;
		f[tid][j].z -= mag * dir.z / dist;
		// }
	}
	#pragma omp single // Суммарная сила будет хранится в первой строке – f[0][i]
	{
		for (int i = 0; i < n; i++) {
			for (int tid = 1; tid < nthreads; tid++) {
				f[0][i].x += f[tid][i].x;
				f[0][i].y += f[tid][i].y;
				f[0][i].z += f[tid][i].z;
			}
		}
	}
}

void move_particles_var1(struct particle *p, struct particle *f, struct particle *v, float *m, int n,
double dt)
{
	for (int i = 0; i < n; i++) {
		struct particle dv = {
			.x = f[i].x / m[i] * dt,
			.y = f[i].y / m[i] * dt,
			.z = f[i].z / m[i] * dt,
		};
		struct particle dp = {
			.x = (v[i].x + dv.x / 2) * dt,
			.y = (v[i].y + dv.y / 2) * dt,
			.z = (v[i].z + dv.z / 2) * dt,
		};

		v[i].x += dv.x;
		v[i].y += dv.y;
		v[i].z += dv.z;
		p[i].x += dp.x;
		p[i].y += dp.y;
		p[i].z += dp.z;
		f[i].x = f[i].y = f[i].z = 0;
	}
}

void move_particles_var2(struct particle *p, struct particle *f, struct particle *v, float *m, int n,
double dt)
{
	#pragma omp for nowait
	for (int i = 0; i < n; i++) {
		struct particle dv = {
			.x = f[i].x / m[i] * dt,
			.y = f[i].y / m[i] * dt,
			.z = f[i].z / m[i] * dt,
		};
		struct particle dp = {
			.x = (v[i].x + dv.x / 2) * dt,
			.y = (v[i].y + dv.y / 2) * dt,
			.z = (v[i].z + dv.z / 2) * dt,
		};

		v[i].x += dv.x;
		v[i].y += dv.y;
		v[i].z += dv.z;
		p[i].x += dp.x;
		p[i].y += dp.y;
		p[i].z += dp.z;
		f[i].x = f[i].y = f[i].z = 0;
	}
}

void move_particles_var5(struct particle *p, struct particle *f[], struct particle *v, float *m, int n,
double dt)
{
	#pragma omp for
	for (int i = 0; i < n; i++) {
		struct particle dv = {
			.x = f[0][i].x / m[i] * dt,
			.y = f[0][i].y / m[i] * dt,
			.z = f[0][i].z / m[i] * dt,
		};

		struct particle dp = {
			.x = (v[i].x + dv.x / 2) * dt,
			.y = (v[i].y + dv.y / 2) * dt,
			.z = (v[i].z + dv.z / 2) * dt,
		};

		v[i].x += dv.x;
		v[i].y += dv.y;
		v[i].z += dv.z;
		p[i].x += dp.x;
		p[i].y += dp.y;
		p[i].z += dp.z;
	}
}

int main()
{
	int x[5] = {1, 2, 4, 6, 8};
	
	for (int i = 0; i < 5; i++) {
		printf("Num threads = %d\n", x[i]);
		double t_total_all = omp_get_wtime();
		for (int j = 0; j < 5; j++)
		{
			double ttotal, tinit = 0, tforces = 0, tmove = 0;

			ttotal = omp_get_wtime();
			// int n = (argc > 1) ? atoi(argv[1]) : 10;
			int n = 50;
			// char *filename = "out_file.txt";
			char *filename = NULL;
			tinit = -omp_get_wtime();

			struct particle *p = malloc(sizeof(*p) * n); // Положение частицы
			struct particle *f = malloc(sizeof(*f) * n); // Сила, действующая на каждую частицу
			struct particle *v = malloc(sizeof(*v) * n); // Скорость частицы
			float *m = malloc(sizeof(*m) * n); // Масса частицы

			for (int i = 0; i < n; i++) {
				p[i].x = rand() / (float)RAND_MAX - 0.5;
				p[i].y = rand() / (float)RAND_MAX - 0.5;
				p[i].z = rand() / (float)RAND_MAX - 0.5;
				v[i].x = rand() / (float)RAND_MAX - 0.5;
				v[i].y = rand() / (float)RAND_MAX - 0.5;
				v[i].z = rand() / (float)RAND_MAX - 0.5;
				m[i] = rand() / (float)RAND_MAX * 10 + 0.01;
				f[i].x = f[i].y = f[i].z = 0;
			}
			tinit += omp_get_wtime();

			double dt = 1e-5;

			#if 0 // var 1
			for (double t = 0; t <= 1; t += dt) { // Цикл по времени (модельному)
				tforces -= omp_get_wtime();
				calculate_forces_var1(p, f, m, n, x[i]); // Вычисление сил – O(N^2)
				tforces += omp_get_wtime();
				tmove -= omp_get_wtime();
				move_particles_var1(p, f, v, m, n, dt); // Перемещение тел O(N)
				tmove += omp_get_wtime();
			}
			#endif

			#if 0 // var 2
			#pragma omp parallel num_threads(x[i])
			{
				for (double t = 0; t <= 1; t += dt) { // Цикл по времени (модельному)
					#pragma omp atomic
					tforces -= omp_get_wtime();
					calculate_forces_var2(p, f, m, n); // Вычисление сил – O(N^2)
					#pragma omp atomic
					tforces += omp_get_wtime();
					#pragma omp atomic
					tmove -= omp_get_wtime();
					#pragma omp barrier
					move_particles_var2(p, f, v, m, n, dt); // Перемещение тел O(N)
					#pragma omp atomic
					tmove += omp_get_wtime();
					#pragma omp barrier
				}
			}
			#endif

			#if 0 // var 3
			locks = malloc(sizeof(omp_lock_t) * n);
			for (int i = 0; i < n; i++) {
				omp_init_lock(&locks[i]);
			}

			#pragma omp parallel num_threads(x[i])
			{
				for (double t = 0; t <= 1; t += dt) { // Цикл по времени (модельному)
					#pragma omp atomic
					tforces -= omp_get_wtime();
					calculate_forces_var3(p, f, m, n); // Вычисление сил – O(N^2)
					#pragma omp atomic
					tforces += omp_get_wtime();
					#pragma omp atomic
					tmove -= omp_get_wtime();
					#pragma omp barrier
					move_particles_var2(p, f, v, m, n, dt); // Перемещение тел O(N)
					#pragma omp atomic
					tmove += omp_get_wtime();
					#pragma omp barrier
				}
			}

			#endif

			#if 0 //var 4

			#pragma omp parallel num_threads(x[i])
			{
				for (double t = 0; t <= 1; t += dt) { // Цикл по времени (модельному)
					#pragma omp atomic
					tforces -= omp_get_wtime();
					calculate_forces_var4(p, f, m, n); // Вычисление сил – O(N^2)
					#pragma omp atomic
					tforces += omp_get_wtime();
					#pragma omp atomic
					tmove -= omp_get_wtime();
					#pragma omp barrier
					move_particles_var2(p, f, v, m, n, dt); // Перемещение тел O(N)
					#pragma omp atomic
					tmove += omp_get_wtime();
					#pragma omp barrier
				}
			}

			#endif

			#if 1 //var 5

			struct particle **f_ = malloc(sizeof(struct particle) * 8);
			// printf("omp max threads = %d\n", x[i]);
			for (int i = 0; i < 8; i++) {
				f_[i] = malloc(sizeof(struct particle) * n);
			}

			#pragma omp parallel num_threads(x[i])
			{
				for (double t = 0; t <= 1; t += dt) { // Цикл по времени (модельному)
					tforces -= omp_get_wtime();
					calculate_forces_var5(p, f_, m, n); // Вычисление сил – O(N^2)
					tforces += omp_get_wtime();
					tmove -= omp_get_wtime();
					#pragma omp barrier
					move_particles_var5(p, f_, v, m, n, dt); // Перемещение тел O(N)
					tmove += omp_get_wtime();
					#pragma omp barrier
				}
			}

			#endif

			ttotal = omp_get_wtime() - ttotal;
			printf("# NBody (n=%d)\n", n);
			printf("# Elapsed time (sec): ttotal %.6f\ntinit %.6f\ntforces %.6f\ntmove %.6f\n",
			ttotal, tinit, tforces, tmove);

			if (filename) {
				FILE *fout = fopen(filename, "w");
				if (!fout) {
					fprintf(stderr, "Can't save file\n");
					exit(EXIT_FAILURE);
				}
				for (int i = 0; i < n; i++) {
					fprintf(fout, "%.3f %.3f %.3f\n", p[i].x, p[i].y, p[i].z);
				}
				fprintf(fout, "\n");
				fclose(fout);
			}

			// free(m);
			// free(v);
			// for (int i = 0; i < x[i]; i++)
			// 	free(f_[i]);
			// free(f_);
			// free(f);
			// free(p);
			// free(locks);
			printf("\n");

		}
		t_total_all = (omp_get_wtime() - t_total_all) / 5;
		printf("\t\t\t\tTOTAL = %.5lf\n", t_total_all);
	}

	return 0;
}