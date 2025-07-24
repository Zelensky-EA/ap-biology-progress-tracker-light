
import React, { useState } from 'react';
import { courseData } from './constants';
import { useStudentProgress } from './hooks/useStudentProgress';
import UnitView from './components/UnitView';
import CourseProgress from './components/CourseProgress';
import VocabularyModal from './components/VocabularyModal';
import UnitSelection from './components/UnitSelection';

const App: React.FC = () => {
  const [selectedUnitId, setSelectedUnitId] = useState<string | null>(null);
  const [selectedVocabTerm, setSelectedVocabTerm] = useState<string | null>(null);
  
  const { 
    progress, 
    updateProgress, 
    calculateICanProgress, 
    unitProgressCalculations, 
    overallCourseProgress 
  } = useStudentProgress();

  const selectedUnit = courseData.units.find(u => u.id === selectedUnitId);

  if (!selectedUnit) {
    return (
      <main className="min-h-screen flex items-center justify-center bg-slate-50 text-slate-800">
        <UnitSelection 
          units={courseData.units} 
          onSelectUnit={setSelectedUnitId}
          overallProgress={overallCourseProgress}
          unitProgress={unitProgressCalculations} 
        />
      </main>
    );
  }

  return (
    <div className="min-h-screen flex flex-col md:flex-row">
      <VocabularyModal term={selectedVocabTerm} onClose={() => setSelectedVocabTerm(null)} />
      
      <aside className="w-full md:w-64 bg-slate-100/80 backdrop-blur-lg border-r border-slate-200 p-4 md:fixed md:h-full overflow-y-auto">
        <div className="mb-4">
          <button onClick={() => setSelectedUnitId(null)} className="w-full text-left group">
            <h1 className="text-2xl font-bold text-slate-900 group-hover:text-cyan-600 transition-colors">AP Biology</h1>
          </button>
        </div>
        <nav>
          <ul>
            {courseData.units.map(unit => (
              <li key={unit.id}>
                <button
                  onClick={() => setSelectedUnitId(unit.id)}
                  className={`w-full text-left p-3 my-1 rounded-md transition-colors text-sm font-medium ${
                    selectedUnitId === unit.id
                      ? 'text-slate-900'
                      : 'text-slate-500 hover:bg-slate-200/50 hover:text-slate-900'
                  }`}
                  style={{
                    backgroundColor: selectedUnitId === unit.id ? unit.color + '22' : undefined,
                    borderLeft: selectedUnitId === unit.id ? `4px solid ${unit.color}` : '4px solid transparent',
                  }}
                >
                  {unit.name}
                </button>
              </li>
            ))}
          </ul>
        </nav>
      </aside>

      <main className="flex-1 md:ml-64">
        <UnitView
            unit={selectedUnit}
            studentProgress={progress}
            onStatusChange={updateProgress}
            calculateICanProgress={calculateICanProgress}
            onVocabClick={(term) => setSelectedVocabTerm(term)}
            unitProgress={unitProgressCalculations[selectedUnit.id] || 0}
        />
      </main>
    </div>
  );
};

export default App;