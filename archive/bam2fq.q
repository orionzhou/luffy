import org.broadinstitute.sting.queue.QScript

class MyScript extends QScript {
  @Input(INPUT, OUTPUT_DIR, OUTPUT_PER_RG=false)
  var scriptInput: File = _
  def script = {
    var myCL = new MyCommandLine
    myCL.myInput = scriptInput // Example variable input
    myCL.myOutput = new File("/path/to/output") // Example hardcoded output
    add(myCL)
  }
}
