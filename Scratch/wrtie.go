// Writing files in Go follows similar patterns to the
// ones we saw earlier for reading.

package main

import (
    "io/ioutil"
    "strconv"
)

func check(e error) {
    if e != nil {
        panic(e)
    }
}

func main() {

    // To start, here's how to dump a string (or just
    // bytes) into a file.
    x, y, z, R, alpha := 25.75, 33.25, 48.5, 36.199, 180.0
    d1 := []byte("title\nScoops3D example A (Arai and Tagyo 1985 example1) homogeneous material properties\nlengthunits   ceeunits  gammaunits\nm   kPa   kN/m^3\nwater\nno\nnmat\n1\nlnum   cee   phi   gamt\n1   41.65   15   18.82\neq\n0\nmethod\nB\nsrch\nsingle\nxcen ycen zcen rad angle\n"+ strconv.FormatFloat(x,'E',-1, 64)+" " + strconv.FormatFloat(y,'E',-1, 64) +" "+ strconv.FormatFloat(z,'E',-1, 64)+" "+ strconv.FormatFloat(R,'E',-1, 64) +" "+ strconv.FormatFloat(alpha,'E',-1, 64) +"\nremove   foscut\nM   10\nisqout\n0\nirelfos\n0\nicritlattice\n0\nisubsurf zfrac\n0   1\nDEM file\nAraiTagyo/input/emb20DEM.asc\noutput directory\nOutput/")
    err := ioutil.WriteFile("test.scp", d1, 0644)
    check(err)
}