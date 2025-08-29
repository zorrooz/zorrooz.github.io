<template>
  <div class="container py-5">
    <div class="row justify-content-center">
      <div class="col-lg-10 col-xl-8">

        <!-- 页面标题 -->
        <div class="text-center mb-6">
          <h1 class="display-5 fw-bold text-primary">关于我</h1>
        </div>

        <!-- 自我介绍 -->
        <section class="mb-5">
          <h2 class="h4 fw-semibold text-dark border-bottom pb-2 mb-3">自我介绍</h2>
          <p class="text-secondary mb-0">
            这是自我介绍的占位符...
          </p>
        </section>

        <!-- 教育经历 -->
        <section class="mb-5">
          <h2 class="h4 fw-semibold text-dark border-bottom pb-2 mb-3">教育经历</h2>
          <div class="ms-2">
            <div v-for="edu in education" :key="edu.degree" class="mb-4">
              <h3 class="h6 fw-bold text-dark mb-1">{{ edu.degree }}</h3>
              <p class="text-muted mb-1 small">
                <strong>{{ edu.institution }}</strong> · {{ edu.duration }}
              </p>
              <p class="text-secondary small mb-0">{{ edu.desc }}</p>
            </div>
          </div>
        </section>

        <!-- 专业技能 -->
        <section class="mb-5">
          <h2 class="h4 fw-semibold text-dark border-bottom pb-2 mb-3">专业技能</h2>
          <div class="ms-2">
            <div v-for="skill in skills" :key="skill.name" class="mb-3">
              <h3 class="h6 fw-bold text-dark mb-1">{{ skill.name }}</h3>
              <p class="text-muted mb-1 small">{{ skill.level }}</p>
              <p class="text-secondary small mb-0">{{ skill.tools }}</p>
            </div>
          </div>
        </section>

        <!-- 项目经历 -->
        <section class="mb-5">
          <h2 class="h4 fw-semibold text-dark border-bottom pb-2 mb-3">项目经历</h2>
          <div class="ms-2">
            <div v-for="project in projects" :key="project.title" class="mb-4">
              <h3 class="h6 fw-bold text-dark mb-1">{{ project.title }}</h3>
              <p class="text-muted mb-1 small">
                <strong>{{ project.role }}</strong> · {{ project.duration }}
              </p>
              <p class="text-secondary small mb-0">{{ project.desc }}</p>
            </div>
          </div>
        </section>

        <!-- 荣誉与奖项 -->
        <section class="mb-5">
          <h2 class="h4 fw-semibold text-dark border-bottom pb-2 mb-3">荣誉与奖项</h2>
          <div class="ms-2">
            <div v-for="(award, index) in awards" :key="index" class="mb-3">
              <h3 class="h6 fw-bold text-dark mb-0">{{ award.title }}</h3>
              <p class="text-muted mb-0 small">
                {{ award.org }} · {{ award.year }}
              </p>
            </div>
          </div>
        </section>

        <!-- 联系我 -->
        <section class="mb-5">
          <h2 class="h4 fw-semibold text-dark border-bottom pb-2 mb-3">联系我</h2>
          <div class="text-center">
            <div class="d-flex flex-wrap justify-content-center align-items-center gap-3 text-muted small">
              <template v-for="contact in contacts" :key="contact.label">
                <a :href="contact.link" target="_blank" rel="noopener"
                  class="d-flex align-items-center text-decoration-none text-secondary hover-primary">
                  <i :class="contact.icon" class="me-1"></i>
                  <span>{{ contact.value }}</span>
                </a>
              </template>
            </div>
          </div>
        </section>

        <!-- 结尾标语 -->
        <div class="text-center mt-5 pt-4 border-top">
          <p class="text-secondary fs-6 mb-0">
            <i class="fas fa-star me-1"></i>
            感谢您的关注！
          </p>
        </div>

      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: 'AboutMeView',
  data() {
    return {
      education: [
        {
          degree: '本科 · 生物信息学',
          institution: '兰州大学',
          duration: '2021 – 2025',
          desc: 'GPA: 3.8/5.0'
        }
      ],
      skills: [
        {
          name: '编程语言',
          level: '熟练掌握',
          tools: 'Python, R, Bash'
        },
        {
          name: '结构生物学工具',
          level: '熟练',
          tools: 'PyMOL, AlphaFold2, ChimeraX'
        },
        {
          name: '数据分析与可视化',
          level: '熟练',
          tools: 'Pandas, NumPy, Matplotlib, ggplot2'
        }
      ],
      projects: [
        {
          title: '这是一个项目',
          role: '项目角色',
          duration: '2024.03 – 2024.06',
          desc: '项目具体描述...'
        }
      ],
      awards: [
        {
          title: '优秀毕业生',
          org: '兰州大学',
          year: '2025'
        }
      ],
      contacts: [
        {
          label: '邮箱',
          value: 'zorrooz@163.com',
          link: 'mailto:zorrooz@163.com',
          icon: 'fas fa-envelope'
        },
        {
          label: 'GitHub',
          value: 'zorrooz',
          link: 'https://github.com/zorrooz',
          icon: 'fab fa-github'
        }
      ]
    }
  }
}
</script>

<style scoped>
.hover-primary:hover {
  color: var(--bs-primary) !important;
}

a.d-flex:hover i {
  transform: scale(1.2);
  transition: transform 0.2s ease;
}
</style>
