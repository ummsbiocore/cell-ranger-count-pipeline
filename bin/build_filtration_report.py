#!/usr/bin/env python3

import sys, argparse, subprocess
from textwrap import dedent

def main(args):

	build_report(args)

def build_report(args):
	output = []
	output.append(print_header(args))
	output.append(setup())
	output.append(read_data())
	output.append(overall_summary())
	output.append(by_criteria())
	output.append(session_info())
	write_output(output)
	run_markdown()

def print_header(args):

	return(dedent(
	'''
	---
	title: "Filtration Summary Report"
	author: "Via Scientific"
	output: 
	  html_document:
	    toc: true
	    toc_float:
	      toc_collapsed: true
	      toc_depth: 3
	    code_folding: hide
	params:
	  input_directory: "{}"
	---
	''').format(args.input_directory).strip())

def setup():
	
	return(dedent(
	'''
	```{r setup, include=FALSE}
	library(dplyr)
	library(ggplot2)
	library(tidyr)
	library(DT)
	```'''))

def read_data():
	
	return(dedent(
	'''
	```{r read_data}
	input_files = list.files(path = params$input_directory, pattern = "filter_summary.tsv$", ignore.case = TRUE)
	
	data = data.frame(signature = factor(),	count=numeric(),	sample=factor())
	for (file in input_files) {
	  data = data %>% add_row(read.delim(paste0(params$input_directory, '/', file), colClasses = c('factor', 'numeric', 'factor')))
	}
	```'''))

def overall_summary():
	
	return(dedent(
	'''
	# Overall {.tabset}
	
	## By Count
	
	```{r overall_count}
	colors = c("Pass"="#6697c2", "Filtered Out"="#fe7f65")
	
	overall = data %>% 
	  mutate(status = ifelse(signature == '1111111', "Pass", "Filtered Out")) %>%
	  group_by(status, sample) %>%
	  summarise(count = sum(count), .groups='drop') %>%
	  group_by(sample) %>%
	  mutate(percent = count / sum(count) * 100)
	    
	ggplot(overall, aes(x=sample, y=count, fill=status, label=count)) +
	  theme_classic() +
	  theme(legend.title = element_blank(),
	        axis.title.x = element_blank(),
	        axis.ticks.x = element_blank(),
	        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) +
	  labs(y = "Number of Cells") +
	  scale_fill_manual(values=colors) +
	  geom_bar(stat='identity', position='stack') +
	  geom_text(position= position_stack())
	```
	
	## By Percentage
	
	```{r overall_percentage}
	ggplot(overall, aes(x=sample, y=percent, fill=status, label=round(percent, 1))) +
	  theme_classic() +
	  theme(legend.title = element_blank(),
	        axis.title.x = element_blank(),
	        axis.ticks.x = element_blank(),
	        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) +
	  labs(y = "Percentage of Cells") +
	  scale_fill_manual(values=colors) +
	  geom_bar( position='stack', stat='identity') +
	  geom_text(position = position_stack())
	```
	
	## Table
	
	```{r}
	overall_summary = overall %>% pivot_wider(names_from = 'status', values_from = count, id_cols = !percent) %>% select(Sample=sample, Pass, `Filtered Out`)
	datatable(overall_summary,
	  rownames = FALSE,
	  extensions = 'Buttons',
	  options=list(
	    dom = 'lftBipr',
	    buttons = list(
	      list(extend = 'csvHtml5', text='Download', filename = "filtration_results", extension='.tsv', fieldBoundary='', fieldSeparator='\\t')
	    )
	  )
	)
	write.table(overall_summary, file='overall_filtration_summary.tsv', quote=FALSE, row.names = FALSE, sep='\\t')
	```'''))

def by_criteria():
	
	return(dedent(
	'''
	# By Criteria
	
	Sample and Count column represent the number of cells that fit into each group. For the critera columns:
	
	1 indcates passing the respective critera
	
	0 indicates failing the respective criteria
	
	```{r by_criteria}
	by_category = data %>% 
	  separate(signature, c("extra", "Min Gene", "Max Gene", "Min UMI", "Max UMI", "Mitochondria", "Ribosome", "Doublet"), "") %>%
	  select(Sample=sample, Count=count, `Min Gene`, `Max Gene`, `Min UMI`, `Max UMI`, Mitochondria, Ribosome, Doublet) %>%
	  arrange(-Count)
	
	sketch = htmltools::withTags(table(
	  class = 'display',
	  thead(
	      th("Sample"),
	      th(style = "border-right: solid 2px;", "Count"),
	      th("Min Gene"),
	      th("Max Gene"),
	      th("Min UMI"),
	      th("Max UMI"),
	      th("% Mitochondrial"),
	      th("% Ribosomal"),
	      th("Doublet Detection"),
	    )
	  )
	)
	
	datatable(by_category,
	  rownames = FALSE,
	  extensions = 'Buttons',
	  options=list(
	    dom = 'lftBipr',
	    buttons = list(
	      list(extend = 'csvHtml5', text='Download', filename = "filtration_criteria_results", extension='.tsv', fieldBoundary='', fieldSeparator='\\t')
	    ),
	  container=sketch
	  )
	) %>% formatStyle(c(2), `border-right` = "solid 2px")
	
	write.table(by_category, file='by_criteria_summary.tsv', quote=FALSE, row.names = FALSE, sep='\\t')
	```'''))

def session_info():

	return(dedent(
	'''
	```{r, session_info, echo=FALSE, results='asis'}
	cat('# Session Info {.tabset .tabset-pills} \\n')
	cat('## Hide\\n')
	cat('## Show\\n')
	sessionInfo()
	```'''))

def write_output(output):

	with open('filtration_summary_report.Rmd', 'w') as out:
		out.write('\n'.join(output)) 

def run_markdown():

	cmd = "Rscript -e 'rmarkdown::render(\"filtration_summary_report.Rmd\", \"html_document\")'"
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	print(out.decode())
	print(err.decode(), file=sys.stderr)

def parseArguments():
	parser = argparse.ArgumentParser(prog="build_QC_report.py", description='', usage='%(prog)s [options]')
	input_args = parser.add_argument_group('Input')
	input_args.add_argument('-i', '--input-directory', required=True, help='Input directory', metavar='', dest='input_directory')
	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)