teste <- function() {
	tabs <- c("Servidores", "Coletores", "Sensores", "Parametros", "Outros")
	dadosConfig <- readConfig("Config Casseta.xls", tabs)
	names(dadosConfig) <- tabs
	
	dadosLer <- prepConfig(dadosConfig)
	dadosLidos <- lerDados(dadosLer)
	}
	
readConfig <- function(file,tabs) {
	require(xlsx)
	lapply(1:5, (function(i) { read.xlsx(file, tabs[i]) } ) )
	}
	
prepConfig <- function(cfg) {
	dadosServidores = cfg[[1]][c(1, 3)]
	dadosColetores  = cfg[[2]][c(1, 3)]
	dadosSensores   = cfg[[3]][c(1, 5, 6)]
	dadosParametros = cfg[[4]]
	dadosOutros     = cfg[[5]]
	aux1 <- merge(dadosSensores, dadosColetores, 
			      by = intersect("Coletor", "Id"))
	aux2 <- merge(aux1, dadosServidores, 
		          by = intersect("Servidor", "Id"))
	names(aux2)[1] <- "Id"
	lerSensores <- aux2[c(1, 2, 3, 5, 7)]
	list(lerSensores, dadosParametros, dadosOutros)	
	}
	
lerDados <- function(cfg) {
	lidoSensores <- lerSensores(cfg[[1]])
	lidoParametros <- lerParametros(lidoSensores, cfg[[2]])
	lidoOutros <- lerOutros(cfg[[3]])
	list(lidoSensores, lidoParametros, lidoOutros)
	}
	
lerSensores <- function(cfg) {
	servTab <- split( cfg, "Servidor")
	
	}