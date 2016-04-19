package main

import (
	"fmt"
	"math/rand"
	"math"
	"scoops3d/Part"
	"time"
)

func main() {
	trial := 1000                //Timeframe
	// FOS := make([]float64, trial)
	// StDev := make([]float64, trial)
	counter := 0.0
	counter2 := 0.0

	// Arai & Tagyo
	// x, y, z, R, alpha :=  26.2264, 32.4998, 46.855, 34.4848, 179.82
	alpha := 179.82

	// Mount St.Helens
	// x, y, z, R, alpha:= 561733.0, 5.11809e+06, 3783.04, 1654.15, 126.583
	

	suksd, phiksd, sukmean, phikmean, gamma, su11, phi11 := Part.SoilP()
	
	RRI := 100.0
	
	t := time.Now()
	fmt.Println(t)

	for x := 25.0; x < 30.0; x=x+1 {
		for y := 20.0; y < 40.0; y=y+1 {
			for z := 30.0; z < 50.0; z=z+1{
				for R := 20.0; R < 40.0; R=R+1{
					// fmt.Println(x, y, z, R, alpha)
					FOS := make([]float64, 0)
					for i := 0; i < trial; i++{			
						A := rand.NormFloat64() 
						B := rand.NormFloat64()
						su11[0] = sukmean[0] + suksd[0]*A
						phi11[0] = phikmean[0] + phiksd[0]*B

						// fmt.Println(i)

						// score, err := Part.Partmain(x, y, z, R, alpha, su11, phi11, gamma)
						score, err := Part.Partmain(x, y, z, R, alpha, 7, su11, phi11, gamma)

						if err == false && score < 100.0{
							FOS = append(FOS, score)
							// FOS[i] = score
							counter = counter + 1.0

							if score < 1.0{
								counter2 = counter2 + 1.0
							}
						}
					}

					mean := Part.SUM(FOS)/float64(len(FOS))
					
					StDev := make([]float64, len(FOS))

					for i := 0; i < len(FOS); i++ {
						StDev[i] = (FOS[i]-mean)*(FOS[i]-mean)
					}

					stdev := math.Sqrt(Part.SUM(StDev)/float64(len(StDev)))
					RI := (mean - 1) / stdev
					
					if RI < RRI {
						RRI = RI + 0.0
						fmt.Println(x, y, z, R, alpha)
						fmt.Println("MEAN:", mean, "St.Dev:", stdev)
						fmt.Println("RI:", RRI)
						fmt.Println("PF:", (counter2/float64(len(FOS)))*100.0)
					}
				}
			}
		}
	}
	tt := time.Now()
	fmt.Println(tt)

	fmt.Println("//////////////////////////FINISH MPI 7//////////////////////////////////////")
}


