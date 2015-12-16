package main

import (
    "bufio"
    "fmt"
    "log"
    "strings"
    "os"
    "strconv"
)

func main() {
    file, err := os.Open("test_out.txt")
    if err != nil {
        log.Fatal(err)
    }
    defer file.Close()

    important:= "1st"
    important2:= "2nd"
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
    fmt.Println(important2)
    important3 := strings.Fields(important2)
    fmt.Println(important3[5])
    FOS, _:=(strconv.ParseFloat(important3[5], 64))
    fmt.Println(FOS)

    if err := scanner.Err(); err != nil {
        log.Fatal(err)
    }
}