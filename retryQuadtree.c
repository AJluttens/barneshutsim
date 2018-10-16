#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

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

// QuadTreeNode structure
typedef struct QuadTreeNode
{   
    // pointer to array of particle structures
    struct particle *particles;

    // we'll need to know how many particles are in a quadrangle
    int numParticles;
    
    // width of the current quadrangle
    double width;

    // four children, northeast - northwest - southwest - southeast
    struct QuadTreeNode *ne;
    struct QuadTreeNode *nw;
    struct QuadTreeNode *sw;
    struct QuadTreeNode *se;
} QuadTreeNode;

// QuadTree structure - A tree node pointer "root" containing the address of a QuadTreeNode structure
struct QuadTree
{
    struct QuadTreeNode *root;
};

// function to create a new QuadTree, point its root to NULL
struct QuadTree* quad_tree_init()
{
    struct QuadTree *qt = (struct QuadTree*)malloc(sizeof (struct QuadTree));
    qt->root = NULL;
    return qt;
}

// function to create a new node, its children point to NULL.
// array of particles structures are parsed and given to the particles attribute
struct QuadTreeNode* quad_tree_node_init(struct particle Particles[])
{
    // a pointer "newNode"
    struct QuadTreeNode *newNode = (struct QuadTreeNode*)malloc(sizeof (struct QuadTreeNode));
    if (newNode)
    {   
        newNode->particles = Particles;
   
        // make new quadrangles point to NULL -> capping of QuadTree
        newNode->ne = NULL;
        newNode->nw = NULL;
        newNode->sw = NULL;
        newNode->se = NULL;
    }
    return newNode;
}

// free a node of the quadtree
void quad_tree_node_free(struct QuadTreeNode *node)
{
    if (node)
    {
        quad_tree_node_free(node->ne);
        quad_tree_node_free(node->nw);
        quad_tree_node_free(node->sw);
        quad_tree_node_free(node->se);
        free(node);
    }
}

// free the entire quadtree
void quad_tree_free(struct QuadTree *qt)
{
    // recursive strategy
    if (qt)
    {
        quad_tree_node_free(qt->root);
        free(qt);
    }
}

// get the total mass of a quadrangle
double CalculateTotalMass(struct particle particles[], int N)
{
    int u;
    double Mass = 0.0;
    for(u=0;u<N;u++)
    {
        Mass += particles[u].mass;
    }
    return Mass;
}

typedef struct centerOfMass{
    double rx;
    double ry;
} com_t;

// get the coordinates of the center of mass of a quadrangle
com_t CalculateCenter(struct centerOfMass center, struct particle particles[], double totalMass, int N)
{
    double x = 0.0, y = 0.0;
    int p;

    for(p=0;p<N;p++)
    {
        x += particles[p].rx*particles[p].mass;
        y += particles[p].ry*particles[p].mass;
    }
    x = x/totalMass;
    y = y/totalMass;
    center.rx = x;
    center.ry = y;
    return center;
}

//
typedef struct Euclidean{
    double dr;
    double dx;
    double dy;
} distance_t;

// function to calculate the Euclidean distance between a center of mass and a particle
distance_t CalculateDistance(struct Euclidean distance, struct centerOfMass center, struct particle particle){
    
    distance.dx = (center.rx-particle.rx);
    distance.dy = (center.ry-particle.ry);
    distance.dr = sqrt(distance.dx*distance.dx+distance.dy*distance.dy);
    return distance;
}

// function to check whether to traverse further down the quadtree or not 
int CheckThetaStatus(double threshold, double width, double distance){
    if (width/(double)distance > threshold)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// FIRST ISSUE IS SITUATED HERE, WHEN DO WE STOP?
// QuadTree subdivider function, recursive strategy
// Give the function a starting node, which will be the root of the tree
void quad_tree_node_subdivide(struct QuadTreeNode *node, int N, struct particle particles[], double x1, double x2, double y1, double y2)
{
    if (node)
    {   
        // initialize arrays for quadrangles containing at most N particles
        // sadly static
        struct particle se_system[N];
        struct particle sw_system[N];
        struct particle ne_system[N];
        struct particle nw_system[N];

        // size of the current quadrangle
        double dx, dy;
        double width = (x2-x1);
        
        dx = 0.5*width;
        dy = 0.5*(y2-y1);

        // initialize quadrangle counters and counter 's' to loop over N particles
        int s, p_se_count=0, p_sw_count=0, p_ne_count=0, p_nw_count=0;

        // initialize coordinates for particles
        // we use the particle's coordinates to check in which quadrangle it should be placed
        double px, py;

        // loop over particles in the parent quadrangle
        for(s=0;s<N;s++){
            // retrieve coordinates for particle 's'
            px = particles[s].rx;
            py = particles[s].ry;

            // check in which quadrangle the particle is situated
            if(px >= x1 && px <= x1+dx && py >= y1 && py <= y1+dy){
                sw_system[p_sw_count] = particles[s];
                p_sw_count += 1;
                //printf("SOUTHWEST %d\n", p_sw_count);
                continue; // make sure particle 's' is only placed in one quadrangle
                }

            if(px >= x1+dx && px <= x2 && py >= y1 && py <= y1+dy){
                se_system[p_se_count] = particles[s];
                p_se_count += 1;
                //printf("SOUTHEAST %d\n", p_se_count);
                continue; // make sure particle 's' is only placed in one quadrangle
                }

            if(px >= x1 && px <= x1+dx && py >= y1+dy && py <= y2){
                nw_system[p_nw_count] = particles[s];
                p_nw_count += 1;
                //printf("NORTHWEST %d\n", p_nw_count);
                continue; // make sure particle 's' is only placed in one quadrangle
                }

            if(px >= x1+dx && px <= x2 && py >= y1+dy && py <= y2){
                ne_system[p_ne_count] = particles[s];
                p_ne_count += 1;
                //printf("NORTHEAST %d\n", p_ne_count);
                continue; // make sure particle 's' is only placed in one quadrangle
                }
        }
        
        //printf("SE:%d\tSW:%d\tNE:%d\tNW:%d\n", p_se_count, p_sw_count, p_ne_count, p_nw_count);

        ///////////////////
        //               //                                               
        //   SOUTHEAST   //
        //               //                                            
        ///////////////////
        if(p_se_count >= 2){ /* More than one particle in this quadrangle -> subdivide again */
            node->se = quad_tree_node_init(se_system);
            memcpy(node->se->particles, se_system, p_se_count);
            
            node->se->numParticles = p_se_count;

            node->se->width = width;
            /* recursive strategy */
            quad_tree_node_subdivide(node->se, p_se_count, se_system, x1+dx, x2, y1, y1+dy);
        }

        ///////////////////
        //               //                                               
        //   SOUTHWEST   //
        //               //                                            
        ///////////////////
        if(p_sw_count >= 2){ /* More than one particle in this quadrangle -> subdivide again */
            node->sw = quad_tree_node_init(sw_system);
            memcpy(node->sw->particles, sw_system, p_sw_count);

            node->sw->numParticles = p_sw_count;

            node->sw->width = width;
            /* recursive strategy */
            quad_tree_node_subdivide(node->sw, p_sw_count, sw_system, x1, x1+dx, y1, y1+dy);
        }

        ///////////////////
        //               //                                               
        //   NORTHEAST   //
        //               //                                            
        ///////////////////
        if(p_ne_count >= 2){ /* More than one particle in this quadrangle -> subdivide again */
            node->ne = quad_tree_node_init(ne_system);
            memcpy(node->ne->particles, ne_system, p_ne_count); 
            
            node->ne->numParticles = p_ne_count;

            node->ne->width = width;
            /* recursive strategy */
            quad_tree_node_subdivide(node->ne, p_ne_count, ne_system, x1+dx, x2, y1+dy, y2);
        }
        
        ///////////////////
        //               //                                               
        //   NORTHWEST   //
        //               //                                            
        ///////////////////
        if(p_nw_count >= 2){ /* More than one particle in this quadrangle -> subdivide again */
            node->nw = quad_tree_node_init(nw_system);
            node->nw->particles = nw_system; 
            
            node->nw->numParticles = p_nw_count;

            node->nw->width = width;
            /* recursive strategy */
            quad_tree_node_subdivide(node->nw, p_nw_count, nw_system, x1, x1+dx, y1+dy, y2);
        }
    }
}

void update(struct QuadTree *qt, struct particle particles[], int N, double x1, double x2, double y1, double y2)
{
    quad_tree_node_free(qt->root);
    qt->root = quad_tree_node_init(particles);
    quad_tree_node_subdivide(qt->root, N, particles, x1, x2, y1, y2);
}

typedef struct accelerations{
    double x;
    double y;
} acc_t;

/* Function to get the count of leaf nodes in a binary tree*/
unsigned int getLeafCount(struct QuadTreeNode* node) 
{ 
  if(node == NULL)        
    return 0; 
  if(node->se == NULL && node->sw==NULL && node->nw==NULL && node->ne==NULL)       
    return 1;             
  else 
    return getLeafCount(node->se)+ 
           getLeafCount(node->sw)+
           getLeafCount(node->ne)+
           getLeafCount(node->nw);       
}

// function to calculate acc_x and acc_y on a particle
acc_t CalculateAccelerations(struct particle particle, struct QuadTreeNode *node, double thresholdTheta, struct accelerations acc){
    if (node && node->numParticles > 0)
    {
        struct centerOfMass nodeCenter;
        double nodeMass;
        struct Euclidean distance;
        nodeMass = CalculateTotalMass(node->particles, node->numParticles);
        nodeCenter = CalculateCenter(nodeCenter, node->particles, nodeMass, node->numParticles);
        distance = CalculateDistance(distance, nodeCenter, particle);

        if (CheckThetaStatus(thresholdTheta, node->width, distance.dr) == 1)
        {
            /* recursive strategy */
            CalculateAccelerations(particle, node->se, thresholdTheta, acc);
            CalculateAccelerations(particle, node->sw, thresholdTheta, acc);
            CalculateAccelerations(particle, node->ne, thresholdTheta, acc);
            CalculateAccelerations(particle, node->nw, thresholdTheta, acc);
        }
        else
        {
            double acceleration;
            // Plummer sphere models
            acceleration = nodeMass/cube(distance.dr+0.001);
            acc.x += acceleration*distance.dx;
            acc.y += acceleration*distance.dy;       
        }
    }
    else
    {
        printf("end of tree!\n");
    }
    return acc;
}



// driver function
int main(int argc, char const *argv[])
{

    ///////////////////////////////////////////////////////////////
    //                                                           //
    //              Reading file and data section                //
    //                                                           //
    ///////////////////////////////////////////////////////////////

    // throw error message if invalid numbers of arguments
    if (argc != 7){
        printf("Invalid number of arguments arguments.\nUsage: ./galsim N filename nsteps delta_t theta graphics\n");
        return 1;
        }

    // assign the arguments to variables
    const int N = atoi(argv[1]);
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const double theta = atof(argv[5]);
    const int graphics = atoi(argv[6]);

    // assign constants once, take epsilon as 0.001 in equation, then it is no longer a variable to call
    const double G = 100/N;

    // you assign the length of the array containing particle structures
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

    ///////////////////////////////////////////////////////////////
    //                                                           //
    //             QuadTree initialization section               //
    //                                                           //
    ///////////////////////////////////////////////////////////////

    // create a new empty QuadTree
    struct QuadTree *qt = quad_tree_init();
    // pass array of particles to the root of the tree
    qt->root = quad_tree_node_init(oldsystem);
    qt->root->numParticles = N;
    qt->root->width = 1.0;
    // subdivide the root recursively until every quadrangle contains only one particle
    // initial rectangular box with width 1
    quad_tree_node_subdivide(qt->root, N, oldsystem, 0.0, 1.0, 0.0, 1.0);

    // get new system ready
    struct particle newsystem[N];
    memcpy (newsystem, oldsystem, sizeof (oldsystem));

    int u, h;
    double init_rx, init_ry, init_mass, init_vx, init_vy;

    ///////////////////////////////////////////////////////////////
    //                                                           //
    //                     Simulation section                    //
    //                                                           //
    ///////////////////////////////////////////////////////////////

    for(h=0;h<nsteps;h++){

        printf("%d\n", getLeafCount(qt->root));

        for(u=0;u<N;u++)
        {
            // get initial coordinates and velocities
            init_rx = oldsystem[u].rx;
            init_ry = oldsystem[u].ry;
            init_vx = oldsystem[u].vx;
            init_vy = oldsystem[u].vy;

            struct accelerations acc;
            acc.x = acc.y = 0.0;
            acc = CalculateAccelerations(qt->root->particles[u], qt->root, theta, acc);

            acc.x *= -G;
            acc.y *= -G;

            // get the new velocities
            init_vx += delta_t*acc.x;
            init_vy += delta_t*acc.y;
            
            // get the new positions
            init_rx += delta_t*init_vx;
            init_ry += delta_t*init_vy;

            // store the new positions and velocities in the particle oldsystem structure
            newsystem[u].rx = init_rx;
            newsystem[u].ry = init_ry;
            newsystem[u].vx = init_vx;
            newsystem[u].vy = init_vy;
        }

        // only copy coordinates and momenta
        memcpy (oldsystem, newsystem, sizeof (oldsystem));

        // update the tree
        update(qt, oldsystem, N, 0.0, 1.0, 0.0, 1.0);

    }

    // free allocated memory for quadtree
    quad_tree_free(qt);

    ///////////////////////////////////////////////////////////////
    //                                                           //
    //                  Writing output section                   //
    //                                                           //
    ///////////////////////////////////////////////////////////////

    FILE *outfile = fopen("assignment4.result.gal", "wb");
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