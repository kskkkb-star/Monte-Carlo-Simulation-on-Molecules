#include<bits/stdc++.h>
#include<time.h>
using namespace std;
#define ll long long
#define newl "\n"

ll molecules = 700;
double box_dim = 10.0; // Length = Breadth = Height = 10 --> simulation space dimensions
ll moves = 200000; // Total number of iterations to do (counting the accepted moves only) 
double kT = 1;
double sigma = 1; 
double epsilon = 1;
double cutoff = 3*sigma;

//random number between left and right --> left always = 0, right = 700(max.)
double random_number(double left, double right)
{
    double range = right-left;
    return (range)*(double(rand())/RAND_MAX) + left;  
    //RAND_MAX is the maximum random number that can be generated = 2147483647 
    //rand() will generate a random number
    //therefore double(rand())/RAND_MAX will generate a random value from (0-1)
    //returning a random number which will be multiple of the range as required
}

//minimum image convention
double min_img(double x)
{
    while(x<-5.0)
        x += 10.0;

    while(x>5.0)
        x -= 10.0;

    return x;
}

/*

The minimum image convention is a technique used in molecular simulations to reduce computational 
overhead by considering only the nearest periodic image of each molecule. This is because in periodic simulations,
the simulation box is replicated infinitely in all directions, creating periodic images of each atom. 
By using the minimum image convention,the simulation can effectively reduce the number of molecules 
that need to be considered for each interaction.

*/

//periodic boundary condition
double pbc(double x)
{
    while(x>10.0)
        x -= 10.0;

    while(x<0)
        x += 10.0;

    return x;
}

/* 

To simulate a bulk system in a finite box by imposing periodicity along the edges of the box. 
This means that when a particle crosses one face of the simulation box, it reappears on the opposite face, 
as if the box were replicated infinitely in all three dimensions. This creates an infinite, 
periodic array of boxes that represent the bulk system. 
PBC can be implemented using different boundary conditions, such as the minimum image convention

#DIFFERENCE BETWEEN MIC AND PBC

In summary, PBC is a technique to simulate a bulk system in a finite simulation box, 
while MIC is a technique to handle interactions between particles in a periodic system, 
by considering only the nearest periodic image of each particle.

*/

//initial energy calculation
double init_energy_calc(double sim_box[][3])
{
    double energy = 0;
    ll i,j;
    for(i=0;i<molecules;i++)
    {
        for(j=i+1;j<molecules;j++)
        {
            double delx = sim_box[j][0] - sim_box[i][0]; 
            double dely = sim_box[j][1] - sim_box[i][1];
            double delz = sim_box[j][2] - sim_box[i][2];

            delx = min_img(delx);  
            dely = min_img(dely);
            delz = min_img(delz);
            /*
            ensuring that if the point j wrt point i is displaced more than -5 to +5; 
            the image of the point j is taken instead
            */ 

            double r_sq = (delx*delx) + (dely*dely) + (delz*delz); 
            //square of distance between two points = (atoms/molecules)

            if(r_sq>cutoff*cutoff) // only calculating for molecules within the cutoff radius
                continue;

            double sigma_sq = sigma*sigma;
            double sr_sq = (sigma_sq/r_sq);
            double a = pow(sr_sq,3);
            double b = a*a;
            double LJ = 4.0*epsilon*(b-a); //U(r) = 4ε[(σ/r)^12 - (σ/r)^6]
            energy += LJ;
            /*
            calculating the LJ Potential Energy between the two points
            (from point 1 with point 2,3....n and then point 2 with point 1,3,4...n and so on)
            */ 
        }
    }
    return energy;
}

//energy after each displacement
double energy_after_disp(double sim_box[][3],double old_x, double old_y, double old_z, double prev_energy,ll random_molecule) 
// random_atom is the generated random_atom which was displaced and its new interactions with the previously fixed moecules gives the new LJ potential energy
{
    double old_interactions = 0.0;
    double new_interactions = 0.0;
    ll i;
    for(i=0;i<molecules;i++)
    {
        if(i!=random_molecule)
        {
            double delx = sim_box[random_molecule][0] - sim_box[i][0];
            double dely = sim_box[random_molecule][1] - sim_box[i][1];
            double delz = sim_box[random_molecule][2] - sim_box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx*delx) + (dely*dely) + (delz*delz);

            if(r_sq>cutoff*cutoff)
                continue;

            double sigma_sq = sigma*sigma;
            double sr_sq = (sigma_sq/r_sq);
            double a = pow(sr_sq,3);
            double b = a*a;
            double LJ = 4.0*epsilon*(b-a); //U(r) = 4ε[(σ/r)^12 - (σ/r)^6]
            new_interactions += LJ;
        }
    }

    for(i=0;i<molecules;i++)
    {
        if(i!=random_molecule)
          {
            double delx = old_x - sim_box[i][0];
            double dely = old_y - sim_box[i][1];
            double delz = old_z - sim_box[i][2];

            delx = min_img(delx);
            dely = min_img(dely);
            delz = min_img(delz);

            double r_sq = (delx*delx) + (dely*dely) + (delz*delz);

            if(r_sq>cutoff*cutoff)
                continue;

            double sigma_sq = sigma*sigma;
            double sr_sq = (sigma_sq/r_sq);
            double a = pow(sr_sq,3);
            double b = a*a;
            double LJ = 4.0*epsilon*(b-a); //U(r) = 4ε[(σ/r)^12 - (σ/r)^6]
            old_interactions += LJ;
        }
    }

    return prev_energy - old_interactions + new_interactions;
}

int main() 
{
    ios::sync_with_stdio(0);
    cin.tie(0);
    srand(time(0));
    double sim_box[molecules][3]; 
    // a 2D array storing information of the points with atoms = atom number and array idicies 0=x,1=y,2=z
    ll i,j;
    vector<double>energy;
    vector<double> index {0,0,0};
    
    
    //starting with a fixed configuration
    for(i=0;i<molecules;i++)
    {
        for(j=0;j<3;j++)
        {
            sim_box[i][j] = (int32_t)((index[j]+0.5)*(10/9)); 
            /*
            using the factor of 10/9 to evenly space out the distribution of the atoms 
            in the box of L=10*10*10 and 700 atoms
            */ 

        }
        index[0] = index[0]+1; // moving x-coordinate x to x+1 as we move from i to i+1 atom
        if(index[0] == 9) //shifting to y-coordinate space and moving y to y+1 as we move from i to i+1 atom
        {
          index[0] = 0;
          index[1] = index[1]+1;
          if(index[1]==9)// repeating the similar process for z to z+1
          {
            index[1]=0;
            index[2] = index[2]+1;
          }
        }
    }

    //pushing the initial finite probability configuration
    energy.push_back(init_energy_calc(sim_box)); 
    cout<<"Move: 0 Energy: "<<energy.back()<<newl; 
    
    ll accept = 0;    
    
    //begining with the simulation
    ll total_attempts = 0; // New counter for total moves attempted
    while(1)
    {
        //selecting a random atom
        total_attempts ++; // New counter for total moves attempted
        ll random_molecule = (ll)(random_number(0,molecules));

        //storing old coordinates for using, if rejected
        double old_x = sim_box[random_molecule][0];
        double old_y = sim_box[random_molecule][1];
        double old_z = sim_box[random_molecule][2];

        //giving random displacement
        sim_box[random_molecule][0] += random_number(0,1.0) - 0.5; // xn = x0 + dleta*(ranf - 0.5)
        sim_box[random_molecule][1] += random_number(0,1.0) - 0.5;
        sim_box[random_molecule][2] += random_number(0,1.0) - 0.5;

        //applying pbc to new coordinates since the random displacement can go beyond the box dimensions
        sim_box[random_molecule][0] = pbc(sim_box[random_molecule][0]);
        sim_box[random_molecule][1] = pbc(sim_box[random_molecule][1]);
        sim_box[random_molecule][2] = pbc(sim_box[random_molecule][2]);

        double new_energy = energy_after_disp(sim_box,old_x,old_y,old_z,energy.back(),random_molecule);

        //energy.back() -> last valid configuration's energy
        double energy_change = new_energy - energy.back();

        //energy change less than 0 -> finite probability
        //so we accept the move
        if(new_energy<=energy.back())
        {
            energy.push_back(new_energy);
            
            accept++;
            //if(accept%5==0)
                cout<<"Move: "<<accept<<" Energy: "<<energy.back()<<newl;
        }
        else
        {
            double prob = exp(-energy_change/kT);

            //calling a random number between 0 and 1
            double ran = random_number(0.0,1.0);

            //if random number <= probability term -> finite probability
            //accept the move 
            //To try not to fall into a local minima and stop the loop only when we obtain the global minima
            if(ran<=prob)
            {
                energy.push_back(new_energy);
                
                accept++;
                //if(accept%5==0)
                    cout<<"Move: "<<accept<<" Energy: "<<energy.back()<<newl;
            }

            //else reject the move and restore the old configuration
            else
            {
                sim_box[random_molecule][0] = old_x;
                sim_box[random_molecule][1] = old_y;
                sim_box[random_molecule][2] = old_z;
            }
        }
        if(accept==moves)
        break;
    }

    double acceptance_rate = (double)accept / total_attempts;
    cout << "Acceptance Rate: " << acceptance_rate * 100 << "%" << newl;
    return 0;
}
