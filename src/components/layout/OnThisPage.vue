<!-- OnThisPage.vue -->
<template>
  <!-- 容器：可选固定高度用于滚动，或自然高度 -->
  <div style="max-height: 300px; overflow-y-auto" class="flex-shrink-0">
    <!-- 标题 -->
    <div class="fw-bold fs-5 text-dark mb-3">On This Page</div>

    <!-- 目录列表 -->
    <ul class="list-unstyled mb-0">
      <li v-for="(section, index) in toc" :key="index" class="mb-2">
        <!-- 二级：章节标题 -->
        <div class="text-dark">{{ section.title }}</div>

        <!-- 三级：子小节，缩进 + 灰色 -->
        <ul v-if="section.children" class="list-unstyled ms-3 mb-0 text-secondary">
          <li v-for="(sub, sIndex) in section.children" :key="sIndex" class="small-text">
            {{ sub.title }}
          </li>
        </ul>
      </li>
    </ul>
  </div>
</template>

<script>
export default {
  name: 'OnThisPage',
  data() {
    return {
      toc: [
        {
          title: "Introduction",
          children: [
            { title: "Getting Started" },
            { title: "Core Concepts" },
            { title: "Setup Guide" }
          ]
        },
        {
          title: "Main Features",
          children: [
            { title: "Real-time Updates" },
            { title: "Data Binding" },
            { title: "Component Reuse" }
          ]
        },
        {
          title: "Advanced Topics",
          children: [
            { title: "State Management" },
            { title: "Performance Tips" },
            { title: "Testing Strategies" }
          ]
        },
        {
          title: "API Reference",
          children: [
            { title: "Props & Events" },
            { title: "Slots" },
            { title: "Lifecycle Hooks" }
          ]
        },
        {
          title: "Conclusion",
          children: [
            { title: "Best Practices" },
            { title: "Next Steps" }
          ]
        }
      ]
    };
  }
};
</script>

<style scoped>
/* 确保滚动条可见（可选美化） */
div[style*="overflow-y-auto"] {
  scrollbar-width: thin;
}

div[style*="overflow-y-auto"]::-webkit-scrollbar {
  width: 6px;
}

div[style*="overflow-y-auto"]::-webkit-scrollbar-thumb {
  background-color: #adb5bd;
  border-radius: 3px;
}

div[style*="overflow-y-auto"]::-webkit-scrollbar-track {
  background: transparent;
}

/* 三级标题字体大小与二级一致（移除 .small） */
.small-text {
  font-size: 1rem;
  /* 与 body 一致，Bootstrap 默认大小 */
  line-height: 1.5;
}
</style>
