import sst
import sst.macro

opticalLatency = "0ps"
smallLatency = "1ms"

def makeUniLink(linkType,srcComp,srcId,srcPort,dstComp,dstId,dstPort,outLat=None,inLat=None):
	if not outLat : outLat = inLat
	if not inLat: inLat = outLat
	if not outLat: sys.exit("must specify at least one latency for link")

	linkName = "%s%d:%d->%d:%d" % (linkType,srcId,srcPort,dstId,dstPort)
	link = sst.Link(linkName)
	portName = "output %d %d" % (srcPort, dstPort)
	srcComp.addLink(link,portName,outLat)
	portName = "input %d %d" % (srcPort, dstPort)
	dstComp.addLink(link,portName,inLat)

def makeUniNetworkLink(srcComp,srcId,srcPort,dstComp,dstId,dstPort,outLat=None,inLat=None):
	makeUniLink("network",srcComp,srcId,srcPort,dstComp,dstId,dstPort,outLat,inLat)

def makeBiLink(linkType,comp1,id1,port1,comp2,id2,port2,outLat=None,inLat=None):
	makeUniLink(linkType,comp1,id1,port1,comp2,id2,port2,outLat,inLat)
	makeUniLink(linkType,comp2,id2,port2,comp1,id1,port1,outLat,inLat)

class Interconnect:
	def __init__(self, params):
		topologiesWithOptics = [
			"flexfly",
			"flexfly_topology"
		]
		self.params = params
		self.simplifiedSwitch = None
		self.system = sst.macro.System(params)
		self.num_nodes = self.system.numNodes()
		self.num_switches = self.system.numSwitches()
		self.switches = [0]*self.num_switches
		self.nodes = [0]*self.num_nodes #just expands the size of an initial array with 1 element 0 to size 1 * num_nodes
		topologyName = (self.params["topology"])["name"]
		if topologyName in topologiesWithOptics:
			self.containsOptics = True
			self.optical_topology_params = self.params["topology"]
			self.num_groups = int(self.optical_topology_params["groups"])
			self.switches_per_group = int(self.optical_topology_params["switches_per_group"])
			self.num_elec_switches = self.num_groups * self.switches_per_group
			self.nodes_per_switch = int(self.optical_topology_params["nodes_per_switch"])
			self.optical_switch_radix = int(self.optical_topology_params["optical_switch_radix"])
			self.optical_switch_radix = self.num_groups
			if self.optical_topology_params.has_key("link_multiplicity"):
				self.link_multiplicity = int(self.optical_topology_params["link_multiplicity"])
			else:
				self.link_multiplicity = 1
			# really just equals to h in Kim's literature (num of global links per switch)
			self.num_nodes = self.num_groups * self.switches_per_group * self.nodes_per_switch
			self.optical_switches_per_group = self.switches_per_group * self.link_multiplicity / self.optical_switch_radix
			if self.optical_switches_per_group % self.optical_switch_radix > 0:
				self.optical_switch_radix += 1
		else:
			self.containsOptics = False

	def buildEndpoints(self):
		nodeParams = self.params["node"]
		modelName = nodeParams["model"]
		for i in range(self.num_nodes): 
			self.nodes[i] = sst.Component("Node %d" % i, "macro.%s" % (modelName + "_node"))
			self.nodes[i].addParams(macroToCoreParams(nodeParams))
			self.nodes[i].addParam("id", i)       
		return

	def buildOpticalSwitches(self):
		switchParams = self.params["switch"]
		switchName = "flexfly_optical" + "_switch"
		#switchName = switchParams["model"] + "_switch"
		for index in range(self.switches_per_group):
			optical_switch_id = index + self.num_elec_switches
			switch = sst.Component("Switch %d" % optical_switch_id, "macro.%s" % switchName)
			switch.addParam("id" , optical_switch_id)
			switch.addParams(macroToCoreParams(switchParams))
			switch.addParam("switch_type" , "optical")
			switch.addParam("optical_switch_radix", self.opticalSwitchRadix)
			switch.addParam("num_electrical_switches", self.num_elec_switches)
			self.switches[optical_switch_id] = switch
		return

	def buildElectricalSwitches(self):
		switchParams = self.params["switch"]
		switchName = "flexfly_electrical" + "_switch"
		#totalSwitchesPerGroup = self.optical_switches_per_group + self.switches_per_group
		for i in range(self.num_elec_switches):
			#if self.containsOptics and i % (totalSwitchesPerGroup) >= self.switches_per_group:
				#continue
			switch = sst.Component("Switch %d" % i, "macro.%s" % switchName)
			switch.addParam("id" , i) 
			switch.addParams(macroToCoreParams(switchParams))
			switch.addParam("switch_type" , "electrical")
			switch.addParam("total_radix", self.switches_per_group + self.nodes_per_switch)
			switch.addParam("switches_per_group", self.switches_per_group)
			switch.addParam("num_groups", self.num_groups);
			self.switches[i] = switch 
		return

	def buildSimplifiedTopology(self, latency):
		self.simplifiedSwitch = sst.Component("Switch 0" , "macro.my_logp_switch")
		self.simplifiedSwitch.addParam("id" , 90) 
		self.simplifiedSwitch.addParam("optical_bandwidth" , "10Gb/s") 
		self.simplifiedSwitch.addParam("electrical_bandwidth" , "1Mb/s") 
		#switch.addParam("")
		for i in range(self.num_nodes):
			makeBiLink("network", self.simplifiedSwitch, 90, i, self.nodes[i], i, 0, latency, latency)
		return

	def latency(self, params):
		if params.has_key("latency"):
			return params["latency"]
		elif params.has_key("send_latency"):
			return params["send_latency"]
		else:
			sys.exit("need link latency in parameters")

	def buildTopology(self):
		switchParams = self.params["switch"]
		for i in range(self.num_switches):
			linkParams = switchParams["link"]
			connections = self.system.switchConnections(i)
			srcSwitch = self.switches[i]
			lat = self.latency(linkParams)
			for srcId, dstId, srcOutport, dstInport in connections:
				#print "srcId: %d and dstId: %d srcOutport: %d, dstInport: %d" % (srcId, dstId, srcOutport, dstInport)
				#print self.switches
				print "connecting - srcID: %d -> dstID: %d" % (srcId, dstId)
				dstSwitch = self.switches[dstId]
				linkName = "network %d:%d->%d:%d" % (srcId, srcOutport ,dstId, dstInport)
        		link = sst.Link(linkName)
        		portName = "output %d %d" % (srcOutport, dstInport)
        		srcSwitch.addLink(link, portName, lat)
        		portName = "input %d %d" % (srcOutport, dstInport)
        		dstSwitch.addLink(link, portName, lat)
				#makeUniNetworkLink(srcSwitch, srcId, srcOutport,
				#					dstSwitch, dstId, dstInport, lat)
		return
	## THIS ONE DOES NOT BUILD IT LIKE LOGP
	def buildTopology2(self):
		switchParams = self.params["switch"]
		for i in range(self.num_switches):
			linkParams = switchParams["link"]
			connections = self.system.switchConnections(i)
			srcSwitch = self.switches[i]
			lat = self.latency(linkParams)
			for srcId, dstId, srcOutport, dstInport in connections:
				#print "srcId: %d and dstId: %d srcOutport: %d, dstInport: %d" % (srcId, dstId, srcOutport, dstInport)
				#print self.switches
				print "connecting - srcID: %d port %d -> dstID: %d port %d" % (srcId, srcOutport, dstId, dstInport)
				dstSwitch = self.switches[dstId]
				makeUniLink("network", srcSwitch, srcId, srcOutport, dstSwitch, dstId, dstInport, smallLatency)
				#linkName = "network %d:%d->%d:%d" % (srcId, srcOutport ,dstId, dstInport)
        		#link = sst.Link(linkName)
        		#portName = "output %d %d" % (srcOutport, dstInport)
        		#srcSwitch.addLink(link, portName, lat)
        		#portName = "input %d %d" % (srcOutport, dstInport)
        		#dstSwitch.addLink(link, portName, lat)

		return		

	def buildNodeConnections(self):
		for i in range(self.num_groups * self.switches_per_group):
			switch = self.switches[i]
			for nodeIndex in range(self.nodes_per_switch):
				index = nodeIndex + (i * self.nodes_per_switch)
				node = self.nodes[index]
				switchPortIndex = self.switches_per_group - 1 + 1 + nodeIndex
				makeUniLink("injection",node,index,0,switch,i,switchPortIndex,smallLatency)
				makeUniLink("ejection",switch,i,switchPortIndex,node,index,0,smallLatency)
			#switchId = self.system.nodeToLogPSwitch(i)
			#switchId = i / self.nodes_per_switch
			#linkname = "InjectionLink%d:Node%d->Switch%d" % (linkCnt, i, switchId)
			#link = sst.Link(linkname)
			#portName1 = "srcPort: %d -> dstPort: %d" % (linkCnt,linkCnt)
			#portName2 = "srcPort: %d -> dstPort: %d" % (linkCnt,linkCnt)
			#node = self.nodes[i]
			
			#switch = self.switches[switchId]
			#portName1 = "in-out %d %d" % (sst.macro.NICLogPInjectionPort, sst.macro.SwitchLogPInjectionPort)
			#portName2 = "in-out %d %d" % (i, sst.macro.SwitchLogPInjectionPort)
			#node.addLink(link, portName1, smallLatency)
			#switch.addLink(link, portName2, smallLatency)
			#linkCnt += 1

	def buildNodeConnections2(self):
		for i in range(self.num_groups * self.switches_per_group):
			switch = self.switches[i]
			for nodeIndex in range(self.nodes_per_switch):
				index = nodeIndex + (i * self.nodes_per_switch)
				node = self.nodes[index]
				switchPortIndex = self.switches_per_group + nodeIndex
				linkName = "logP %d -> %d" % (index, i)
				link = sst.Link(linkName)
				print "NICLogInjectionPort: %d and SwitchLogPInjectionPort: %d" % (sst.macro.NICLogPInjectionPort, sst.macro.SwitchLogPInjectionPort)
				portName = "in-out %d %d" % (sst.macro.NICLogPInjectionPort, sst.macro.SwitchLogPInjectionPort)
				#portName = "in-out %d %d" % (0, self.switches_per_group - 1 + nodeIndex)
				node.addLink(link, portName, smallLatency)
				#portName = "in-out %d %d" % (0, self.switches_per_group - 1 + nodeIndex)
				portName = "in-out %d %d" % (i, sst.macro.SwitchLogPInjectionPort)
				switch.addLink(link, portName, smallLatency)


	def buildLogPNetwork(self):
		import re
		nproc = sst.getMPIRankCount() * sst.getThreadCount()
		switchParams = self.params["switch"]
		linkParams = switchParams["link"]
		ejParams = switchParams["ejection"]
		lat = self.latency(ejParams)
		#gotta multiply the lat by 2
		match = re.compile("(\d+[.]?\d*)(.*)").search(lat)
		if not match:
			sys.exit("improperly formatted latency %s" % lat)
		num, units = match.groups()
		num = eval(num) * 2
		lat = "%8.4f%s" % (num,units.strip())
		switches = []
		for i in range(nproc):
			switch = sst.Component("my_logp %d" % i, "macro.my_logp_switch")
			#switch = sst.Component("logp %d" % i, "macro.logp_switch")
			switch.addParams(macroToCoreParams(switchParams))
			switch.addParam("id", i)
			switches.append(switch)
			#NEWLY ADDED BEGIN
			switch.addParam("electrical_bandwidth" , "10Gb/s")
			switch.addParam("optical_bandwidth", "1Gb/s")
			#NEWLY ADDED END
		for i in range(nproc):
			sw_i = switches[i]
			for j in range(nproc):
				sw_j = switches[j]
				if i==j: continue
				linkName = "logPnetwork%d->%d" % (i,j)
				link = sst.Link(linkName)
				portName = "in-out %d %d" % (j, sst.macro.SwitchLogPNetworkPort)
				sw_i.addLink(link, portName, lat)
				portName = "in-out %d %d" % (i, sst.macro.SwitchLogPNetworkPort)
				sw_j.addLink(link, portName, lat)			
		
		for i in range(self.num_nodes):
			injSW = self.system.nodeToLogPSwitch(i)
			ep = self.nodes[i]
			sw = switches[injSW]
			linkName = "logPinjection%d->%d" % (i, injSW)
			link = sst.Link(linkName)
			portName = "in-out %d %d" % (sst.macro.NICLogPInjectionPort, sst.macro.SwitchLogPInjectionPort)
			#portName = "injection %d %d" % (sst.macro.NICLogPInjectionPort, sst.macro.SwitchLogPInjectionPort)
			ep.addLink(link, portName, smallLatency) #put no latency here
			portName = "in-out %d %d" % (i, sst.macro.SwitchLogPInjectionPort)
			#portName = "ejection %d %d" % (i, sst.macro.SwitchLogPInjectionPort)
			sw.addLink(link, portName, smallLatency)

	def makeOneOpticalSwitch(self, switchParams, i):
		opticalSwitchName = "macro." + "flexfly_optical_switch"
		opticalSwitch = sst.Component("Switch %d" % i, opticalSwitchName)
		opticalSwitch.addParam("id" , i) 
		opticalSwitch.addParams(macroToCoreParams(switchParams))
		opticalSwitch.addParam("switch_type" , "optical")
		opticalSwitch.addParam("optical_switch_radix", self.opticalSwitchRadix)
	
	def buildNetworkMonitor(self):
		nodeParams = self.params["node"]
		self.network_manager_node = sst.Component("Node %d" % self*num_nodes, "macro.%s" % ("network_manager" + "_node"))
		self.network_manager_node.addParams(macroToCoreParams(nodeParams))
		id = self.num_nodes
		self.network_manager_node.addParam("id", i);
		return

	def build(self, islogP):
		if not islogP:
			self.buildEndpoints()
			topologyParams = self.params["topology"]
			topologyName = topologyParams["name"]
			self.buildElectricalSwitches()
			if self.containsOptics:
				self.opticalSwitchRadix = topologyParams["optical_switch_radix"]
				self.buildOpticalSwitches()
			#self.buildTopology()
			self.buildTopology2()
			self.buildNodeConnections()
			#self. buildNodeConnections2()
			self.buildLogPNetwork()
		else:
			self.buildEndpoints()
			self.buildSimplifiedTopology(smallLatency)
			self.buildLogPNetwork()


		

## Returns the command line argv in terms of a vector that is 0-indexed
def readCmdLineParams():
	import sys
	return sst.macro.readParams(sys.argv)

def redoSubParams_impl(nsArr, theDict, allParams):
	for key in theDict:
		val = theDict[key]
		if isinstance(val, dict):
			newNsArr = nsArr[:]
			newNsArr.append(key)
			redoSubParams_impl(newNsArr, val, allParams)
		else:
			paramArr = nsArr[:]
			paramArr.append(key)
			newParam = ".".join(paramArr)
			allParams.append((newParam, val))

def macroToCoreParams(theDict):
	allParams = []
	redoSubParams_impl([], theDict, allParams)
	newDict = {}
	for key, val in allParams:
		newDict[key] = val
	return newDict

def setupTopology():
	import sys
	params = readCmdLineParams()
	nodeParams = params["node"]
	swParams = params["switch"]
	
	builtinApps = [
   		"apitest",
       	"global_test",
       	"hello_world",
       	"mpi_coverage",
       	"mpi_ping_all",
       	"mpi_print_nodes",
       	"mpi_topology",
       	"parsedumpi",
       	"sstmac_mpi_testall",
       	"traffic_matrix",
       	"user_app_cxx_empty_main",
       	"user_app_cxx_full_main",
	]

	for appIdx in range(10):
		appKey = "app%d" % (appIdx)
		#nodeParams = params["node"]
		if params.has_key(appKey):
			appParams = params[appKey]
			nodeParams[appKey] = appParams
			appName = appParams["name"]
			if not appName in builtinApps:
				cmd = "import sst.%s" % appName
				exec(cmd)
			del params[appKey]

	debugList = []
	if params.has_key("debug"):
		debugList = params["debug"].strip().split()
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-d" or sys.argv[i] == "--debug":
			debugList.extend(sys.argv[i+1].split(","))

	icParams = {}
	icParams["topology"] = params["topology"]
	nodeParams["interconnect"] = icParams
	if debugList:
		nodeParams["debug"] = " ".join(debugList)
	swParams["topology"] = params["topology"]

	swParams["topology"] = "torus"

	#move every param in the global namespace 
	#into the individal namespaces
	for ns in "node", "switch":
		nsParams = params[ns]
		for key in params:
			val = params[key]
			if isinstance(val, str):
				if not nsParams.has_key(key):
					nsParams[key] = val    
	interconnect = Interconnect(params)
	isLogP = False
	interconnect.build(isLogP)
