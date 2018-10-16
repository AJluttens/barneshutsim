#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// a function to find the cube of a certain number
double cube(double num)
{
    return (num * num * num);
}

// structure for information on particles in the simulation
typedef struct particle {
	double rx;
	double ry;
	double mass;
	double vx;
	double vy;
	double brightness;
    } particle_t;

int main(int argc, char const *argv[])
{	
	// if the input does not match, raise and error and exit
	if (argc != 6){
		printf("Invalid number of arguments arguments.\nUsage: ./galsim N filename nsteps delta_t graphics\n");
		return 1;
		}
	
	// assign the arguments to variables
	const int N = atoi(argv[1]);
	const int nsteps = atoi(argv[3]);
	const double delta_t = atof(argv[4]);
	const int graphics = atoi(argv[5]);

	// assign constants once, take epsilon as 0.001 in equation, then it is no longer a variable to call
	const double G = 100/N;

	// dynamically allocate memory for data from input binary file called 'system'
	struct particle oldsystem[N];

	// read the binary file using a FOR LOOP to discretely assign values to the systems information
	int i,j, body=0;
	double d;
	FILE *infile = fopen(argv[2], "r");
	for(i=0; i<6*N; i = i+6){
		for(j=0;j<6;j++){
			fread(&d, sizeof(double), 1, infile);
			switch(j){
			case 0:
				{
				oldsystem[body].rx = d;
				break;
				}
			case 1:
				{
				oldsystem[body].ry = d;
				break;
				}
			case 2:
				{
				oldsystem[body].mass = d;
				break;
				}
			case 3:
				{
				oldsystem[body].vx = d;
				break;
				}
			case 4:
				{
				oldsystem[body].vy = d;
				break;
				}
			case 5:
				{
				oldsystem[body].brightness = d;
				break;
				}
			}
		}
		body += 1;

	}
	struct particle newsystem[N];
	memcpy (newsystem, oldsystem, sizeof (oldsystem));

	int h, k;
	double init_rx, init_ry, init_mass, init_vx, init_vy, dx, dy, dr, acc, acc_x, acc_y;
	// oldsystem evolving over time
	for(h=0;h<nsteps;h++)
	{
		for(i=0;i<N;i++){
			//O(N^2) algorithm, could have used symmetry here
			init_rx = oldsystem[i].rx;
			init_ry = oldsystem[i].ry;
			init_vx = oldsystem[i].vx;
			init_vy = oldsystem[i].vy;
			init_mass = oldsystem[i].mass;
			for(k=0;k<N;k++){
				// skip self interaction
				if(k!=i){
					dx = init_rx - oldsystem[k].rx;
					dy = init_ry - oldsystem[k].ry;
					// Euclidean distance formula
					dr = sqrt(dx*dx + dy*dy);
					printf("Distance: %lf", dr);					
					// this is not correct -- DONT DO THIS TWICE
					acc = oldsystem[k].mass/cube(dr+0.001);
					acc_x += dx*acc;
					acc_y += dy*acc;
				}
			}
			
			acc_x *= -G;
			acc_y *= -G;
			
			printf("acc: %lf, %lf\n", acc_x, acc_y);	

			// get the new velocities
			init_vx += delta_t*acc_x;
			init_vy += delta_t*acc_y;
			
			// get the new positions
			init_rx += delta_t*init_vx;
			init_ry += delta_t*init_vy;
	
			// store the new positions and velocities in the particle oldsystem structure
			newsystem[i].rx = init_rx;
			newsystem[i].ry = init_ry;
			newsystem[i].vx = init_vx;
			newsystem[i].vy = init_vy;

			// reset the acceleration for the new particle in the loop
			acc_x = 0;
			acc_y = 0;
		}
		// only copy coordinates and momenta
		memcpy (oldsystem, newsystem, sizeof (oldsystem));
	}
	
	FILE *outfile = fopen("result.gal", "wb");
	for(i=0;i<N;i++){
		printf("%lf,%lf\n", oldsystem[i].rx, oldsystem[i].ry);
		fwrite(&oldsystem[i].rx,sizeof(double),1,outfile);
		fwrite(&oldsystem[i].ry,sizeof(double),1,outfile);
		fwrite(&oldsystem[i].mass,sizeof(double),1,outfile);
		fwrite(&oldsystem[i].vx,sizeof(double),1,outfile);
		fwrite(&oldsystem[i].vy,sizeof(double),1,outfile);
		fwrite(&oldsystem[i].brightness,sizeof(double),1,outfile);
	}
	fclose(outfile);

return 0;
}
