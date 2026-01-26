#!/bin/bash
# Generate weekly summary prompt
# Usage: bash scripts/weekly_summary.sh

echo "======================================================================"
echo "WEEKLY SUMMARY TEMPLATE"
echo "======================================================================"
echo ""
echo "Week: $(date '+Week %V, %B %Y')"
echo ""
echo "Copy this template to your work log:"
echo ""
cat << 'TEMPLATE'
### Week [X]: [Date Range]

**Major Accomplishments:**
- 
- 
- 

**Challenges Faced:**
- 
- 

**Lessons Learned:**
- 
- 

**Progress Toward Milestones:**
- [x] Completed item
- [ ] In progress item

**Plans for Next Week:**
- 
- 
- 

**Hours Worked:** X hours

**Reflection:**
[Your thoughts on the week]

TEMPLATE
