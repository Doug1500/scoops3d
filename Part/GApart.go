package Part

import (
	// "fmt"
	"math"
	"io/ioutil"
	"strconv"
	"os/exec"
	"bufio"
	"log"
	"strings"
	"os"
)

func Partmain(x, y, z, R, alpha float64, su, phi, gamma []float64 ) (score float64, err bool){
// func Partmain(x, y, z, R, alpha float64,rank int, su, phi, gamma []float64 ) (score float64, err bool){
// func Partmain(x, y, z, R, alpha, FOS float64, Igood []float64) (score float64, err bool){
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////SINGLE SOIL MATERIAL/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	FileName := "AraiTagyo"
	ASCFile := "emb20DEM.asc"
	material := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)
	// material := "1   41.65   15   18.82"
			 // tag	Cohesion phi gamma	

	// FileName := "StHelens"
	// ASCFile := "sthel_res100mDEM.asc"
	// material := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)
	// material := "1   1000.0   40.0   24.0"
	
	// FileName := "Cone"
	// ASCFile := "cone1000DEM.asc"
	// material := "1   17630   40   21"
	// material := "1   1444   30   20"
	// material := "1   172.9   20   19"

	// FileName := "DonaldGiam"
	// ASCFile := "emb10DEM.asc"
	// material := "1   3   19.6   20"

	d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nno\nnmat\n1\nlnum   cee   phi   gamt\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\noutput directory\nOutput/")

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////MULTIPLE MATERIAL/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	
	// material1 := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)
	// material2 := "2   "+strconv.FormatFloat(su[1],'E',-1, 64)+"   "+strconv.FormatFloat(phi[1],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[1],'E',-1, 64)
	// material3 := "3   "+strconv.FormatFloat(su[2],'E',-1, 64)+"   "+strconv.FormatFloat(phi[2],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[2],'E',-1, 64)

	// material1 := "1   29.4   12   18.82"
	// material2 := "2   9.8   5   18.82"
	// material3 := "3   294   40   18.82" 
	// LayerFile := "emb20layer"

	// FileName := "DonaldGiam"
	// ASCFile := "emb10DEM.asc"
	// LayerFile := "emb10layer"
	// material1 := "1   "+strconv.FormatFloat(su[0],'E',-1, 64)+"   "+strconv.FormatFloat(phi[0],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[0],'E',-1, 64)
	// material2 := "2   "+strconv.FormatFloat(su[1],'E',-1, 64)+"   "+strconv.FormatFloat(phi[1],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[1],'E',-1, 64)
	// material3 := "3   "+strconv.FormatFloat(su[2],'E',-1, 64)+"   "+strconv.FormatFloat(phi[2],'E',-1, 64)+"   "+strconv.FormatFloat(gamma[2],'E',-1, 64)

	// material1 := "1   0   38   19.5"
	// material2 := "2   5.3   23   19.5"
	// material3 := "3   7.2   20   19.5" 
	// LayerFile := "emb10layer"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nno\nnmat\n3\nlnum   cee   phi   gamt\n"+material1+"\n"+material2+"\n"+material3+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   5.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\nlayer file\n"+FileName+"/input/"+LayerFile+"\noutput directory\nOutput/")

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////3D MULTIPLE MATERIAL DATA////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material := "1   -1   -1   -1"
	// MaterialFile := "emb20mat3D.txt"

	// FileName := "DonaldGiam"
	// ASCFile := "emb10DEM.asc"
	// material := "1   -1   -1   -1"
	// MaterialFile := "emb10mat3d.txt"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nno\nstr3d   linterp\n1   0\nnmat\n1\nlnum   cee   phi   gamt\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   5.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\nmaterial properties file\n"+FileName+"/input/"+MaterialFile+"\noutput directory\nOutput/")	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////PIEZOMETER HEAD DATA/////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material := "1   41.65   15   18.82   18.82"
	// PiezoFile := "emb20piezo.asc"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater  gamw\npz   9.81\nnmat\n1\nlnum   cee   phi   gamt\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\npiezometric file\n"+FileName+"/input/"+PiezoFile+"\noutput directory\nOutput/")	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////PRESSURE HEAD DATA///////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FileName := "AraiTagyo"
	// ASCFile := "emb20DEM.asc"
	// material := "1   41.65   15   18.82   18.82"
	// Piezo3DFile := "emb20phead3D.txt"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater  gamw\n3d   9.81\nnmat\n1\nlnum   cee   phi   gamt\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\npressure head file\n"+FileName+"/input/"+Piezo3DFile+"\noutput directory\nOutput/")	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////MATERIAL DATA RU///////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FileName := "Cone"
	// ASCFile := "cone1000DEM.asc"
	// material := "1   1444   30   20   0.2"

	// d1 := []byte("title\nScoops3D example R; Mount Saint Helens\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nru\nnmat\n1\nlnum   cee   phi   gamt   ru\n"+material+"\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10.0\nisqout\n1\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\n"+FileName+"/input/"+ASCFile+"\noutput directory\nOutput/")

	// err1 := ioutil.WriteFile("test"+strconv.Itoa(rank)+".scp", d1, 0644)
	err1 := ioutil.WriteFile("test.scp", d1, 0644)
	Check(err1)

	//exe.go --> 
	exec.Command("/home/yewintun/mygo/src/scoops3d/fortrancode/scoops3d.exe").Run()
	
	//read.go --> read test_out.txt -->FOS
	// file, err3 := os.Open("/home/yewintun/mygo/src/scoops3d/Output/test"+strconv.Itoa(rank)+"_out.txt")
	file, err3 := os.Open("/home/yewintun/mygo/src/scoops3d/Output/test_out.txt")

	if err3 != nil {
		log.Fatal(err3)
	}

	defer file.Close()

	important, important2 := "1st", "2nd"
	counter := 0.0

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		a:=scanner.Text()
		if counter < 2.0{
			if a == "3D POTENTIAL FAILURE"{
				important = a
			}
			if important == "3D POTENTIAL FAILURE"{
				important2 = a
				counter = counter + 1.0
			}
		}
	}
	
	important3 := strings.Fields(important2)
	// fmt.Println(important3)
	if important3[0] == "2nd" {
		score = 1000.0
	}else{
		score, _ =(strconv.ParseFloat(important3[5], 64))
		if score == 0.0{
			score = 1000.0
		}
	}

	if err3 := scanner.Err(); err3 != nil {
		log.Fatal(err3)
	}

	if err3 := scanner.Err(); err3 != nil {
		log.Fatal(err3)
	}

	// if score== 1000.0 || score == 100.0{
	// 	I := []float64{x, y, z, R, alpha}
	// 	score = SetValue(FOS, Igood, I)
	// }

	return score, err
}

func Check(e error) {
	if e != nil {
		panic(e)
	}
}

func SetValue(FOS float64, Igood, I []float64) float64 {
	Ifinal := make([]float64, len(Igood))
	for i := 0; i < len(Igood); i++ {
		Ifinal[i] = (Igood[i]-I[i])*(Igood[i]-I[i])
	}
	score := math.Sqrt(SUM(Ifinal))*1000.0 + FOS
	// fmt.Println(score, I)
	return score
}

func SUM(slice []float64) float64 {
	sum := 0.0
	for i := 0; i < len(slice); i++ {
		sum = slice[i] + sum
	}
	return sum
}

func SoilP() ([]float64, []float64, []float64, []float64, []float64, []float64, []float64) {
	//////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////AriaTagyo////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	// suksd, phiksd := 8.0, 3.0
	suksd := []float64{8.0, 0.0, 0.0}
	phiksd:= []float64{3.0, 0.0, 0.0}

	sukmean := []float64{41.65, 41.65, 41.65}
	phikmean:= []float64{15.0, 15.0, 15.0}
	gamma 	:= []float64{18.82, 18.82, 18.82}
	
	su11 := []float64{41.65, 41.65, 41.65}
	phi11:= []float64{15.0, 15.0, 15.0}

	//////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////Donald&Giam//////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	// suksd, phiksd := 0.000002, 0.000001

	// sukmean := []float64{0.0, 5.3, 7.2}
	// phikmean:= []float64{38.0, 23.0, 20.0}
	// gamma 	:= []float64{19.5, 19.5, 19.5}
	
	// su11 := []float64{0.0, 5.3, 7.2}
	// phi11:= []float64{38.0, 23.0, 20.0}
	
	// //////////////////////////////////////////////////////////////////////////////
	// /////////////////////////////////St.Helens////////////////////////////////////
	// //////////////////////////////////////////////////////////////////////////////
	// sukmean := []float64{1000.0, 1000.0, 1000.0}
	// phikmean:= []float64{40.0, 40.0, 40.0}
	// gamma 	:= []float64{24.0, 24.0, 24.0}

	// suksd, phiksd := 636.3961, 7.071

	// su11 := []float64{1000.0, 1000.0, 1000.0}
	// phi11:= []float64{40.0, 40.0, 40.0}


	return suksd, phiksd, sukmean, phikmean, gamma, su11, phi11
}
