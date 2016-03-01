package main

import (
	"fmt"
	"math/rand"
	"math"
	"scoops3d/Part"
)

func main() {
	trial := 100000                //Timeframe
	// FOS := make([]float64, trial)
	// StDev := make([]float64, trial)
	FOS := make([]float64, 0)
	

	counter := 0.0
	counter2 := 0.0

	// x, y, z, R, alpha :=  24.6543, 34.1117, 50.5277, 37.5, 178.229
	// x, y, z, R, alpha :=  561186.0, 5117810.0, 3590.33, 2001.1, 139.631
	x, y, z, R, alpha := 564126.0, 5.11825e+06, 3837.82, 1970.84, 45.3499

	suksd, phiksd, sukmean, phikmean, gamma, su11, phi11 := Part.SoilP()
	
	for i := 0; i < trial; i++ {
		A := rand.NormFloat64() 
		B := rand.NormFloat64()
		su11[0] = sukmean[0] + suksd*A
		phi11[0] = phikmean[0] + phiksd*B

		fmt.Println(i)

		// score, err := Part.Partmain(x, y, z, R, alpha, su11, phi11, gamma)
		score, err := Part.Partmain(x, y, z, R, alpha, su11, phi11, gamma)

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
	
	fmt.Println(x, y, z, R, alpha)
	fmt.Println("MEAN:", mean, "St.Dev:", stdev)
	fmt.Println("RI:", RI)
	fmt.Println("PF:", (counter2/float64(len(FOS)))*100.0)

	// fmt.Println(FOS)
	// fmt.Println(StDev)

}


