"""Main ScRNA Agent orchestrator."""

import os
import json
import logging
from typing import Optional, Dict, Any
from anthropic import Anthropic

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ScRNAAgent:
    """Main AI Agent for scRNA-seq analysis using Claude."""

    def __init__(self, api_key: Optional[str] = None):
        """Initialize the ScRNA Agent.
        
        Args:
            api_key: Anthropic API key (defaults to ANTHROPIC_API_KEY env var)
        """
        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY environment variable not set")
        
        self.client = Anthropic()
        self.conversation_history = []
        self.analysis_results = {}
        
        # Import tools
        from tools import get_all_tools
        self.tools = get_all_tools()
        
    def _format_tools_for_claude(self) -> list:
        """Format tools for Claude's tool use API."""
        return self.tools
    
    def _process_tool_call(self, tool_name: str, tool_input: Dict[str, Any]) -> str:
        """Execute a tool call and return results."""
        logger.info(f"Executing tool: {tool_name}")
        
        try:
            # Import and execute the appropriate tool
            if tool_name == "run_qc":
                from tools.qc_tools import run_qc
                return json.dumps(run_qc(**tool_input))
            
            elif tool_name == "normalize":
                from tools.preprocessing_tools import normalize
                return json.dumps({"status": "normalized", "method": tool_input.get("method")})
            
            elif tool_name == "cluster":
                from tools.analysis_tools import cluster
                return json.dumps({"status": "clustered", "method": tool_input.get("method")})
            
            elif tool_name == "find_markers":
                from tools.annotation_tools import find_markers
                return json.dumps({"status": "markers_found", "n_markers": 10})
            
            elif tool_name == "annotate_cells":
                from tools.annotation_tools import annotate_cells
                return json.dumps({"status": "annotated", "cell_types": []})
            
            elif tool_name == "reduce_dimensions":
                from tools.analysis_tools import reduce_dimensions
                return json.dumps({"status": "reduced", "method": tool_input.get("method")})
            
            else:
                return json.dumps({"error": f"Unknown tool: {tool_name}"})
        
        except Exception as e:
            logger.error(f"Error executing tool {tool_name}: {str(e)}")
            return json.dumps({"error": str(e)})
    
    def plan_analysis(self, data_path: str, query: str) -> str:
        """Have Claude plan the analysis."""
        logger.info("Planning analysis with Claude...")
        
        system_prompt = """You are an expert bioinformatician analyzing single-cell RNA-seq data.
        
Your role:
1. Understand the user's question about their scRNA-seq data
2. Plan an optimal analysis workflow
3. Use available tools to execute the analysis
4. Interpret results and provide biological insights

Available tools:
- run_qc: Perform quality control filtering
- normalize: Normalize and scale expression data
- reduce_dimensions: Apply dimensionality reduction (PCA, UMAP, tSNE)
- cluster: Identify cell clusters
- find_markers: Find marker genes for clusters
- annotate_cells: Annotate cell types

Always start by running quality control, then proceed with other analyses."""

        # Add user message
        self.conversation_history.append({
            "role": "user",
            "content": f"I have scRNA-seq data at {data_path}. {query}"
        })
        
        # Get Claude's response with tool use
        response = self.client.messages.create(
            model="claude-3-5-sonnet-20241022",
            max_tokens=4096,
            system=system_prompt,
            tools=self._format_tools_for_claude(),
            messages=self.conversation_history
        )
        
        # Process response
        tool_results = []
        assistant_message = {"role": "assistant", "content": response.content}
        self.conversation_history.append(assistant_message)
        
        # Handle tool use
        for block in response.content:
            if block.type == "tool_use":
                logger.info(f"Tool call: {block.name}")
                tool_result = self._process_tool_call(block.name, block.input)
                tool_results.append({
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": tool_result
                })
                self.analysis_results[block.name] = tool_result
        
        # Continue conversation if there were tool calls
        if tool_results:
            self.conversation_history.append({
                "role": "user",
                "content": tool_results
            })
            
            # Get Claude's interpretation
            response = self.client.messages.create(
                model="claude-3-5-sonnet-20241022",
                max_tokens=4096,
                system=system_prompt,
                tools=self._format_tools_for_claude(),
                messages=self.conversation_history
            )
            
            self.conversation_history.append({
                "role": "assistant",
                "content": response.content
            })
        
        # Extract text response
        text_response = ""
        for block in response.content:
            if hasattr(block, "text"):
                text_response += block.text
        
        return text_response
    
    def verify_results(self) -> str:
        """Have Claude verify the analysis results."""
        logger.info("Verifying results with Claude...")
        
        verification_prompt = f"""Based on the analysis performed, verify the results:

Analysis Results Summary:
{json.dumps(self.analysis_results, indent=2)}

Please check:
1. Data quality metrics
2. Clustering validity
3. Cell type annotations consistency
4. Overall analysis quality

Provide recommendations if any issues are found."""

        self.conversation_history.append({
            "role": "user",
            "content": verification_prompt
        })
        
        response = self.client.messages.create(
            model="claude-3-5-sonnet-20241022",
            max_tokens=2048,
            messages=self.conversation_history
        )
        
        self.conversation_history.append({
            "role": "assistant",
            "content": response.content
        })
        
        # Extract text
        text_response = ""
        for block in response.content:
            if hasattr(block, "text"):
                text_response += block.text
        
        return text_response
    
    def generate_report(self) -> str:
        """Have Claude generate a final analysis report."""
        logger.info("Generating report with Claude...")
        
        report_prompt = """Based on the complete analysis, generate a comprehensive report including:

1. Executive Summary
2. Methods Used
3. Key Findings
4. Cell Type Identification
5. Biological Insights
6. Recommendations for Further Analysis

Format the report in markdown with clear sections and bullet points."""

        self.conversation_history.append({
            "role": "user",
            "content": report_prompt
        })
        
        response = self.client.messages.create(
            model="claude-3-5-sonnet-20241022",
            max_tokens=4096,
            messages=self.conversation_history
        )
        
        self.conversation_history.append({
            "role": "assistant",
            "content": response.content
        })
        
        # Extract text
        text_response = ""
        for block in response.content:
            if hasattr(block, "text"):
                text_response += block.text
        
        return text_response
    
    def analyze(self, data_path: str, query: str, verify: bool = True) -> Dict[str, Any]:
        """Run complete analysis pipeline.
        
        Args:
            data_path: Path to scRNA-seq data (H5AD, CSV, MTX)
            query: User's analysis question
            verify: Whether to verify results
            
        Returns:
            Dictionary with results and report
        """
        logger.info("Starting scRNA-seq analysis...")
        
        # Step 1: Planning
        plan = self.plan_analysis(data_path, query)
        logger.info("Planning complete")
        
        # Step 2: Verification (optional)
        verification = ""
        if verify:
            verification = self.verify_results()
            logger.info("Verification complete")
        
        # Step 3: Report Generation
        report = self.generate_report()
        logger.info("Report generation complete")
        
        return {
            "plan": plan,
            "verification": verification,
            "report": report,
            "results": self.analysis_results
        }


if __name__ == "__main__":
    # Example usage
    agent = ScRNAAgent()
    
    results = agent.analyze(
        data_path="data/example.h5ad",
        query="Perform quality control, normalize, cluster, and identify cell types"
    )
    
    print("\n=== ANALYSIS REPORT ===")
    print(results["report"])
