using CalculusWithJulia, Plots,DataFrames, CSV, Statistics
#-------------------------------------------------------------------------------------------------------------------------
# FUNCTION TO CALCULATE MSD FOR ONE TRAJECTORY
function MSDcalculate(x,y,N_Max,N)
    ltrack= N+1
    msd=zeros(N_Max+1)
    for tau in 1:N_Max
        for i in tau+1:ltrack
            msd[tau+1]+=((x[i]-x[i-tau])^2+(y[i]-y[i-tau])^2)/(ltrack-tau)
        end
    end
    return msd
end

#-----------------------------------------------------------------------------------------------------------------
# FUNCTION THAT GENERETAES n TRAJECTORIES, PLOTS FIRST FOUR OF THEM, CALCULATE AVERAGE MSD OF ALL OF THEM
function traj_and_MSD(x0, y0, R::Float64, v::Float64, num_traj::Int64, N, Delta_t::Float64, N_Max)
    graph = plot();
    matrMSD = fill(NaN, N_Max+1, num_traj)

    for i in 1:num_traj
        orientazione = rand()*2*pi
        orientazione = round(orientazione, digits=3)
        abp = initABP( (x0, y0, orientazione), R, v);

        p, t = trajectory( abp, N, Delta_t);

        x = [pi[1] for pi in p]
        y = [pi[2] for pi in p]

        if i <= 4 
            plot!(x,y, range=[-50,50],  title = "ActiveParticle (R=$R µm, v=$v µm/s)", aspect_ratio= :equal, legend=false)
            xlabel!("x [μm]")
            ylabel!("y [μm]")
        end

        matrMSD[1:N_Max+1, i] = MSDcalculate(x,y, N_Max, N)
    end

    return graph, matrMSD 

end


function MSD_single!(x,y,time)
  T= size(time,1)
  R= 2.0
  msd=zeros(T)
  msdv=zeros(T)
    for τ in 1:T
        for i in 1:T-τ
            msd[τ]+= ((x[i+τ]-x[i])^2+(y[i+τ]-y[i])^2)./(T-τ)
            
        end
        msd[τ] = msd[τ]/(R^2)
    end

    k= plot(time/100, msd)
 return k
end

path= "C:\\Users\\j.sharma\\OneDrive - Scuola Superiore Sant'Anna\\P07 Coding\\2023\\12.Dec\\ellipse\\20231221-154400\\R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1\\run1\\"

filename= "20231221-154400 R=2.0 v=10.0 a=50.0 b=25.0 pf=0.1 run1"

pathf= path*filename

f=  pathf*".csv" 

df= CSV.read(f,DataFrame)
#y=df[44:74,:ypos]
#time=df[44:74,:N]
gdf = groupby(df,:N,sort=true)

xgdf= [g[!,:xpos] for g in gdf]
ygdf= [g[!,:ypos] for g in gdf]
tgdf= [g[!,:Time] for g in gdf]
tgdf[1]
T= size(tgdf[1],1)
for i in 1:T
     tt= MSD_single!(xgdf[i],ygdf[i],tgdf[i])
     #@show i
    display(tt)
end


#tt= MSD_single!(x,y,time)