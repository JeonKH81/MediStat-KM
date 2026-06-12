#!/usr/bin/env python3
"""
Kaplan-Meier Survival Analysis Script
NEJM-style visualization with comprehensive statistical output
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from lifelines.utils import median_survival_times
import warnings
warnings.filterwarnings('ignore')

# NEJM-style settings
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Helvetica', 'Liberation Sans'],
    'font.size': 10,
    'axes.linewidth': 1,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'legend.frameon': False,
})

NEJM_COLORS = ['#BC3C29', '#0072B5', '#E18727', '#20854E', '#7876B1', '#6F99AD', '#FFDC91', '#EE4C97']



def load_data(filepath: str) -> pd.DataFrame:
    """Load Excel file."""
    return pd.read_excel(filepath)


def km_analysis(df: pd.DataFrame, time_col: str, event_col: str, group_col: str = None,
                x_interval: int = 30, output_prefix: str = 'km_result',
                group_labels: dict = None) -> dict:
    """
    Perform Kaplan-Meier analysis with NEJM-style output.
    
    Args:
        df: DataFrame with survival data
        time_col: Column name for time-to-event
        event_col: Column name for event (1=event, 0=censored)
        group_col: Column name for grouping variable (optional)
        x_interval: X-axis interval in days (30 or 365)
        output_prefix: Prefix for output files
        group_labels: Dictionary mapping group values to labels (e.g., {0: 'non-DM', 1: 'DM'})
    
    Returns:
        Dictionary with statistical results
    """
    results = {}
    
    # Clean data
    cols = [time_col, event_col] + ([group_col] if group_col else [])
    df_clean = df[cols].dropna()
    
    T = df_clean[time_col]
    E = df_clean[event_col].astype(int)
    
    # Determine x-axis range
    max_time = T.max()
    x_ticks = np.arange(0, max_time + x_interval, x_interval)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    if group_col is None:
        # Single group analysis
        kmf = KaplanMeierFitter()
        kmf.fit(T, E, label='All patients')
        
        # Plot event probability (1 - survival)
        event_prob = 1 - kmf.survival_function_
        ax.step(event_prob.index, event_prob.iloc[:, 0], where='post', color=NEJM_COLORS[0], linewidth=2)
        
        # Confidence interval
        ci_lower = 1 - kmf.confidence_interval_survival_function_.iloc[:, 1]
        ci_upper = 1 - kmf.confidence_interval_survival_function_.iloc[:, 0]
        ax.fill_between(event_prob.index, ci_lower, ci_upper, alpha=0.2, color=NEJM_COLORS[0], step='post')
        
        # Median survival
        median = kmf.median_survival_time_
        ci = median_survival_times(kmf.confidence_interval_)
        results['median_survival'] = {
            'All': {'median': median, 'ci_lower': ci.iloc[0, 0], 'ci_upper': ci.iloc[0, 1]}
        }
        results['n_events'] = {'All': int(E.sum())}
        results['n_total'] = {'All': len(E)}
        
        # Number at risk
        risk_times = x_ticks[x_ticks <= max_time]
        n_at_risk = [kmf.event_table.loc[:t, 'at_risk'].iloc[-1] if t in kmf.event_table.index or len(kmf.event_table.loc[:t]) > 0 else len(T) for t in risk_times]
        results['n_at_risk'] = {'All': dict(zip(risk_times, n_at_risk))}
        
    else:
        # Multiple group analysis
        groups = df_clean[group_col].unique()
        groups = sorted(groups)
        
        # Calculate event rate for each group to determine color assignment
        event_rates = {}
        for group in groups:
            mask = df_clean[group_col] == group
            event_rate = E[mask].sum() / mask.sum()
            event_rates[group] = event_rate
        
        # Sort groups by event rate (highest first gets red color)
        groups_sorted_by_event = sorted(groups, key=lambda x: event_rates[x], reverse=True)
        
        # Create label mapping
        if group_labels is None:
            # Auto-generate labels: "GroupCol=Value" format or use value directly if meaningful
            group_labels = {}
            for g in groups:
                g_val = int(g) if hasattr(g, 'item') else g  # Convert numpy types
                if isinstance(g_val, int) and g_val in [0, 1]:
                    # Binary variable: use descriptive names
                    group_labels[g] = f"non-{group_col}" if g_val == 0 else group_col
                else:
                    group_labels[g] = str(g)
        else:
            # Convert provided labels to match numpy types
            new_labels = {}
            for g in groups:
                g_val = int(g) if hasattr(g, 'item') else g
                if g_val in group_labels:
                    new_labels[g] = group_labels[g_val]
                elif g in group_labels:
                    new_labels[g] = group_labels[g]
                else:
                    new_labels[g] = str(g)
            group_labels = new_labels
        
        # Assign colors: red for highest event rate, blue for lower
        color_map = {}
        for i, group in enumerate(groups_sorted_by_event):
            color_map[group] = NEJM_COLORS[i % len(NEJM_COLORS)]
        
        kmf_dict = {}
        results['median_survival'] = {}
        results['n_events'] = {}
        results['n_total'] = {}
        results['n_at_risk'] = {}
        
        # Plot in order of event rate (highest first) for legend consistency
        for group in groups_sorted_by_event:
            mask = df_clean[group_col] == group
            kmf = KaplanMeierFitter()
            label = group_labels.get(group, str(group))
            kmf.fit(T[mask], E[mask], label=label)
            kmf_dict[group] = kmf
            
            color = color_map[group]
            
            # Plot event probability
            event_prob = 1 - kmf.survival_function_
            ax.step(event_prob.index, event_prob.iloc[:, 0], where='post', 
                   color=color, linewidth=2, label=label)
            
            # Confidence interval
            ci_lower = 1 - kmf.confidence_interval_survival_function_.iloc[:, 1]
            ci_upper = 1 - kmf.confidence_interval_survival_function_.iloc[:, 0]
            ax.fill_between(event_prob.index, ci_lower, ci_upper, alpha=0.15, color=color, step='post')
            
            # Median survival
            median = kmf.median_survival_time_
            ci = median_survival_times(kmf.confidence_interval_)
            results['median_survival'][label] = {
                'median': median, 
                'ci_lower': ci.iloc[0, 0] if not ci.empty else np.nan,
                'ci_upper': ci.iloc[0, 1] if not ci.empty else np.nan
            }
            results['n_events'][label] = int(E[mask].sum())
            results['n_total'][label] = int(mask.sum())
        
        # Log-rank test
        if len(groups) == 2:
            g1, g2 = groups
            lr = logrank_test(T[df_clean[group_col] == g1], T[df_clean[group_col] == g2],
                             E[df_clean[group_col] == g1], E[df_clean[group_col] == g2])
            results['logrank_p'] = lr.p_value
            results['logrank_stat'] = lr.test_statistic
        else:
            lr = multivariate_logrank_test(T, df_clean[group_col], E)
            results['logrank_p'] = lr.p_value
            results['logrank_stat'] = lr.test_statistic
        
        # Unadjusted Hazard Ratio (Cox univariate)
        df_cox = df_clean[[time_col, event_col, group_col]].copy()
        df_cox[group_col] = pd.Categorical(df_cox[group_col]).codes
        
        cph = CoxPHFitter()
        cph.fit(df_cox, duration_col=time_col, event_col=event_col)
        
        hr = np.exp(cph.params_[group_col])
        hr_ci = np.exp(cph.confidence_intervals_.loc[group_col])
        results['unadjusted_hr'] = {
            'HR': hr,
            'CI_lower': hr_ci.iloc[0],
            'CI_upper': hr_ci.iloc[1],
            'p_value': cph.summary.loc[group_col, 'p']
        }
        
        ax.legend(loc='lower right', fontsize=10)
    
    # Formatting
    ax.set_xlabel('Time (days)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Event Probability', fontsize=11, fontweight='bold')
    ax.set_xlim(0, max(x_ticks))
    ax.set_ylim(0, 1)
    ax.set_xticks(x_ticks)
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Number at risk table
    if group_col:
        risk_times = x_ticks[x_ticks <= max_time]
        
        for idx, group in enumerate(groups_sorted_by_event):
            kmf = kmf_dict[group]
            label = group_labels.get(group, str(group))
            n_at_risk_list = []
            for t in risk_times:
                if len(kmf.event_table.loc[:t]) > 0:
                    n_at_risk_list.append(int(kmf.event_table.loc[:t, 'at_risk'].iloc[-1]))
                else:
                    n_at_risk_list.append(int(results['n_total'][label]))
            results['n_at_risk'][label] = dict(zip([int(t) for t in risk_times], n_at_risk_list))
            
            y_pos = -0.12 - idx * 0.05
            ax.text(-0.02, y_pos, label, transform=ax.transAxes, fontsize=9, 
                   fontweight='bold', ha='right', va='center', color=color_map[group])
            for j, (t, n) in enumerate(zip(risk_times, n_at_risk_list)):
                x_pos = t / max(x_ticks)
                ax.text(x_pos, y_pos, str(n), transform=ax.transAxes, fontsize=8, ha='center', va='center')
        
        ax.text(-0.02, -0.07, 'No. at Risk', transform=ax.transAxes, fontsize=9, fontweight='bold', ha='right')
    else:
        # Single-curve number at risk table
        risk_times = x_ticks[x_ticks <= max_time]
        n_at_risk_list = list(results['n_at_risk']['All'].values())
        ax.text(-0.02, -0.07, 'No. at Risk', transform=ax.transAxes, fontsize=9, fontweight='bold', ha='right')
        ax.text(-0.02, -0.12, 'All', transform=ax.transAxes, fontsize=9,
                fontweight='bold', ha='right', va='center', color=NEJM_COLORS[0])
        for t, n in zip(risk_times, n_at_risk_list):
            x_pos = t / max(x_ticks)
            ax.text(x_pos, -0.12, str(int(n)), transform=ax.transAxes, fontsize=8, ha='center', va='center')
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2 if group_col else 0.16)
    
    # Save figure
    fig_path = f'{output_prefix}_km_curve.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    results['figure_path'] = fig_path
    return results




def format_results(results: dict) -> str:
    """Format results as a readable string."""
    output = []
    output.append("=" * 60)
    output.append("KAPLAN-MEIER SURVIVAL ANALYSIS RESULTS")
    output.append("=" * 60)
    
    # Sample size
    output.append("\n[ Sample Size ]")
    for group, n in results['n_total'].items():
        events = results['n_events'][group]
        output.append(f"  {group}: N={n}, Events={events}")
    
    # Median survival
    output.append("\n[ Median Survival Time (95% CI) ]")
    for group, data in results['median_survival'].items():
        median = data['median']
        ci_l = data['ci_lower']
        ci_u = data['ci_upper']
        if np.isnan(median):
            output.append(f"  {group}: Not reached")
        else:
            output.append(f"  {group}: {median:.1f} days ({ci_l:.1f} - {ci_u:.1f})")
    
    # Log-rank test
    if 'logrank_p' in results:
        output.append("\n[ Log-rank Test ]")
        output.append(f"  Chi-square statistic: {results['logrank_stat']:.3f}")
        output.append(f"  P-value: {results['logrank_p']:.4f}")
    
    # Unadjusted HR
    if 'unadjusted_hr' in results:
        hr = results['unadjusted_hr']
        output.append("\n[ Unadjusted Hazard Ratio ]")
        output.append(f"  HR: {hr['HR']:.3f} (95% CI: {hr['CI_lower']:.3f} - {hr['CI_upper']:.3f})")
        output.append(f"  P-value: {hr['p_value']:.4f}")
    
    # Number at risk
    output.append("\n[ Number at Risk ]")
    for group, risk_data in results['n_at_risk'].items():
        times = sorted(risk_data.keys())
        risk_str = ", ".join([f"t={t}: {risk_data[t]}" for t in times[:6]])
        output.append(f"  {group}: {risk_str}...")
    
    output.append("\n" + "=" * 60)
    output.append(f"Figure saved: {results['figure_path']}")
    output.append("=" * 60)
    
    return "\n".join(output)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Kaplan-Meier Survival Analysis')
    parser.add_argument('--file', required=True, help='Path to Excel file')
    parser.add_argument('--time', required=True, help='Time-to-event column name')
    parser.add_argument('--event', required=True, help='Event column name (1=event, 0=censored)')
    parser.add_argument('--group', default=None, help='Group column name (optional)')
    parser.add_argument('--interval', type=int, default=30, help='X-axis interval in days (30 or 365)')
    parser.add_argument('--output', default='km_result', help='Output file prefix')
    parser.add_argument('--labels', nargs='+', default=None, help='Group labels as value:label (e.g., 0:non-DM 1:DM)')

    args = parser.parse_args()

    group_labels = None
    if args.labels:
        group_labels = {}
        for item in args.labels:
            if ':' in item:
                val, label = item.split(':', 1)
                try:
                    val = int(val)
                except ValueError:
                    try:
                        val = float(val)
                    except ValueError:
                        pass
                group_labels[val] = label

    df = load_data(args.file)
    print(f"\nData loaded: {len(df)} rows\n")

    results = km_analysis(df, args.time, args.event, args.group, args.interval, args.output, group_labels)
    print(format_results(results))
    print("\n✅ KM analysis complete.")
