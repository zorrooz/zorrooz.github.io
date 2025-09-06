<!-- Home.vue -->
<template>
  <div class="container py-4">
    <div class="row">

      <!-- 中间主内容区域 -->
      <div class="col-12 col-lg-9 order-1 order-lg-2">
        <div class="container">
          <div class="row">
            <div class="col">
              <!-- <MD :content="markdownContent" /> -->
              <p>abcd 测试</p>
              <h3>文档中心</h3>
              <DocumentList :docs="documentData" :perPage="6" />
            </div>
          </div>
        </div>
      </div>

      <!-- 左侧边栏：移动端显示在最下面 -->
      <div class="col-12 col-lg-3 order-2 order-lg-1">
        <div class="d-flex flex-column align-items-center align-items-lg-end w-100 gap-3">
          <ProfileCard class="w-100" />
          <TagCloud class="w-100" />
        </div>
      </div>

    </div>
  </div>
</template>

<script>
import ProfileCard from '@/components/layout/ProfileCard.vue'
import TagCloud from '@/components/layout/TagCloud.vue'
// import MD from '@/components/layout/MD.vue'
import markdownContent from '@/content/demo.md?raw'
import DocumentList from '@/components/layout/DocumentList.vue'

export default {
  name: 'MainView',
  components: {
    ProfileCard,
    TagCloud,
    // MD,
    DocumentList,
  },
  data() {
    return {
      markdownContent,
      documentData: [
        {
          id: 1,
          title: "Vue 3 入门指南",
          description: "介绍 Vue 3 的基本语法和组合式 API 使用方法。",
          category: ["前端", "Vue", "gblog"],
          date: "2025-03-10",
          tags: ["Vue", "Composition API", "新手"]
        },
        {
          id: 2,
          title: "Pinia 状态管理",
          description: "如何在项目中使用 Pinia 替代 Vuex 进行状态管理。",
          category: "状态管理",
          date: "2025-04-15",
          tags: ["Vue", "Composition API", "新手"]
        },
        {
          id: 3,
          title: "TypeScript 基础语法",
          description: "变量、接口、泛型等核心概念讲解。",
          category: "语言",
          date: "2025-02-20",
          tags: ["Vue", "Composition API", "新手"]
        },
        {
          id: 4,
          title: "RESTful API 设计规范",
          description: "统一接口设计、状态码、资源命名等最佳实践。",
          category: "后端",
          date: "2025-05-01",
          tags: ["Vue", "Composition API", "新手"]
        },
        {
          id: 5,
          title: "Vite 构建优化",
          description: "提升开发服务器启动速度和打包性能。",
          category: "构建工具",
          date: "2025-06-12",
          tags: ["Vue", "Composition API", "新手"]
        },
        {
          id: 6,
          title: "ESLint 与 Prettier 集成",
          description: "代码规范自动化配置方案。",
          category: "工程化",
          date: "2025-07-03",
          tags: ["Vue", "Composition API", "新手"]
        },
        {
          id: 7,
          title: "响应式布局实战",
          description: "使用 Bootstrap 实现多端适配。",
          category: "CSS",
          date: "2025-01-25",
          tags: ["Vue", "Composition API", "新手"]
        },
        {
          id: 8,
          title: "单元测试入门",
          description: "使用 Vitest 编写第一个测试用例。",
          category: "测试",
          date: "2025-08-10",
          tags: ["Vue", "Composition API", "新手"]
        }]
    }
  }
}
</script>

<style scoped></style>
