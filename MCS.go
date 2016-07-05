// /////////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////// CASE STUDY 1 & 3 Mount St Helen and Arai & Tagyo///////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////////////////////
// package main

// import (
// 	"fmt"
// 	"math/rand"
// 	"math"
// 	"scoops3d/Part"
// )

// func main() {
// 	trial := 200                //Timeframe
// 	FOS := make([]float64, 0)

// 	counter := 0.0
// 	counter2 := 0.0

// 	// // Arai & Tagyo
// 	// x, y, z, R, alpha :=  26.22, 32.50, 46.86, 34.48, 179.82

// 	// // MCS + Brute Force
// 	// x, y, z, R, alpha :=  28.0, 32.0, 47.0, 35.0, 180.0

// 	// Mount St.Helens
// 	// x, y, z, R, alpha:= 560847.0, 5117430.0, 3944.1, 2139.68, 149.386
// 	x, y, z, R, alpha := 1044.6, 2824.2, 1404.8, 1518.2, 141.51


// 	suksd, phiksd, sukmean, phikmean, gamma, su11, phi11 := Part.SoilP()
		
// 	for i := 0; i < trial; i++ {
// 		A := rand.NormFloat64() 
// 		B := rand.NormFloat64()
// 		su11[0] = sukmean[0] + suksd[0]*A
// 		phi11[0] = phikmean[0] + phiksd[0]*B

// 		// score, err := Part.Partmain(x, y, z, R, alpha, su11, phi11, gamma)
// 		score, _, err := Part.Partmain(x, y, z, R, alpha, 7, su11, phi11, gamma)
// 		fmt.Println(i, score)

// 		if err == false && score < 100.0{
// 			FOS = append(FOS, score)
// 			// FOS[i] = score
// 			counter = counter + 1.0
			
// 			if score < 1.0{
// 				counter2 = counter2 + 1.0
// 			}
// 		}
// 	}


// 	mean := Part.SUM(FOS)/float64(len(FOS))
	
// 	StDev := make([]float64, len(FOS))

// 	for i := 0; i < len(FOS); i++ {
// 		StDev[i] = (FOS[i]-mean)*(FOS[i]-mean)
// 	}
// 	stdev := math.Sqrt(Part.SUM(StDev)/float64(len(StDev)))
// 	RI := (mean - 1) / stdev
	

// 	fmt.Println(x, y, z, R, alpha)
// 	fmt.Println("MEAN:", mean, "St.Dev:", stdev)
// 	fmt.Println("RI:", RI)
// 	fmt.Println("PF:", (counter2/float64(len(FOS)))*100.0)

// 	// fmt.Println(FOS)
// 	// fmt.Println(StDev)

// }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// CASE STDUY 2 DONALD AND GIAM /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

package main

import (
	"fmt"
	"math/rand"
	"math"
	"scoops3d/Part"
)

func main() {
	trial := 1000                //Timeframe
	// FOS := make([]float64, trial)
	// StDev := make([]float64, trial)
	FOS := make([]float64, 0)
	

	counter := 0.0
	counter2 := 0.0

	// Donald & Giam
	x, y, z, R, alpha :=  15.7945, 52.2399, 54.0327, 23.6607, 167.427

	suksd, phiksd, sukmean, phikmean, gamma, su11, phi11 := Part.SoilP()
	
	
	for i := 0; i < trial; i++ {
		A := rand.NormFloat64() 
		B := rand.NormFloat64()
		C := rand.NormFloat64()
		D := rand.NormFloat64()
		E := rand.NormFloat64()
		

		phi11[2] = phikmean[2] + phiksd[2]*E
		phi11[1] = phikmean[1] + phiksd[1]*D
		phi11[0] = phikmean[0] + phiksd[0]*C
		su11[2] = sukmean[2] + suksd[2]*B
		su11[1] = sukmean[1] + suksd[1]*A

		fmt.Println(i)

		// score, err := Part.Partmain(x, y, z, R, alpha, su11, phi11, gamma)
		score, _ , err := Part.Partmain(x, y, z, R, alpha, 0, su11, phi11, gamma)

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



